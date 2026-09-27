"""Read an ACQ file into a GSData object."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
from astropy import units as un
from astropy.coordinates import EarthLocation, Longitude
from astropy.time import Time
from pygsdata import KNOWN_TELESCOPES, GSData, Telescope
from pygsdata.history import History, Stamp
from pygsdata.readers import gsdata_reader
from pygsdata.select import freq_selector, lst_selector, select_loads, time_selector
from pygsdata.utils import time_concat

from . import _coordinates as _crd
from .read_acq import ACQError, _index_file, _read_spectra, decode_file, encode


def fast_lst_setter(times: Time, loc: EarthLocation):
    """Set the LSTs for the GSData object."""
    years, days, hours, minutes, seconds = np.array(
        [yd.split(":") for yd in times.yday.flatten()]
    ).T

    secs = _crd.tosecs(
        years.astype(int),
        days.astype(int),
        hours.astype(int),
        minutes.astype(int),
        seconds.astype(float),
    )
    gst = _crd.gst(secs) * 12 / np.pi + loc.lon.hour

    return Longitude(gst * un.hour).reshape(times.shape)


LOADS = ("ant", "internal_load", "internal_load_plus_noise_source")
SELECTORS = ("freq_selector", "time_selector", "lst_selector", "load_selector")


@gsdata_reader(select_on_read=True, formats=["acq"])
def read_acq_to_gsdata(
    path: str | Path | Sequence[str | Path],
    telescope: Telescope = KNOWN_TELESCOPES["edges-low"],
    name: str = "{year}-{day}:{hour}:{minute}",
    lst_setter: callable | None = None,
    selectors: dict[str, dict] | None = None,
    **kwargs,
) -> GSData:
    """Read an ACQ file into a GSData object.

    Parameters
    ----------
    path : str or Path
        The path to the file to read. If a sequence of paths is given, the files will be
        concatenated along the time axis.
    telescope_location : str or EarthLocation
        The location of the telescope.
    name : str
        The name of the GSData object. Can include formatting fields for year, day,
        hour, minute and stem (the stem of the input filename).
    lst_setter : callable, optional
        A function that takes the astropy Time object at which the data is defined and
        returns a set of LSTs that will be stored in the GSData object. By default, use
        the astropy function used by pygsdata itself. Set
    selectors : dict, optional
        Selections to apply while reading, so that only the selected spectra (and
        channels) are decoded. Keys may be any of ``freq_selector``,
        ``time_selector``, ``lst_selector`` and ``load_selector``, each mapping to a
        dict of keyword arguments for the corresponding function in
        :mod:`pygsdata.select` (e.g. ``{"lst_selector": {"lst_range": (6, 12)}}``).
        The result is the same as reading everything and then applying
        ``select_freqs``, ``select_times``, ``select_lsts`` and ``select_loads`` in
        that order. When multiple files are given, ``indx`` selections refer to the
        concatenated times. This is also what ``GSData.from_file(...,
        selectors=...)`` uses.
    **kwargs
        Additional keyword arguments to pass to the GSData constructor.
    """
    path = [Path(path)] if isinstance(path, (str, Path)) else [Path(p) for p in path]

    if len(path) == 0:
        raise ValueError("No files given to read")

    if not all(p.exists() for p in path):
        raise FileNotFoundError(f"File {path} does not exist")

    if isinstance(telescope, str):
        try:
            telescope = KNOWN_TELESCOPES[telescope]
        except KeyError as e:
            raise ValueError(
                "telescope must be a Telescope or name of a KNOWN_TELESCOPE, "
                f"got {telescope}. Known 'scopes are {KNOWN_TELESCOPES.keys()}."
            ) from e

    if selectors:
        return _read_selected(
            path,
            telescope=telescope,
            name=name,
            lst_setter=lst_setter,
            selectors=selectors,
            **kwargs,
        )

    # Read the _first_ file to get the metadata
    pants = []
    ploads = []
    plnss = []
    times = []
    ancs = []
    for pth in path:
        _, (pant, pload, plns), anc = decode_file(pth)
        _times = Time(anc.data.pop("times"), format="yday", scale="utc")

        if pant.size == 0:
            continue  # ignore this file, it's empty.

        pants.append(pant)
        ploads.append(pload)
        plnss.append(plns)
        times.append(_times)
        ancs.append(anc)

    if len(pants) == 0:
        raise ACQError(f"No data in any files: {path}")

    pant = np.concatenate(pants, axis=1)
    pload = np.concatenate(ploads, axis=1)
    plns = np.concatenate(plnss, axis=1)
    times = time_concat(times)
    anc = ancs[0]
    anc.data = {k: np.concatenate([a.data[k] for a in ancs]) for k in anc.data}

    return _build_gsdata(
        data=np.array([pant.T, pload.T, plns.T])[:, np.newaxis],
        times=times,
        freqs=anc.frequencies * un.MHz,
        auxiliary_measurements=anc.data,
        path=path,
        telescope=telescope,
        name=_format_name(name, times, path),
        lst_setter=lst_setter,
        **kwargs,
    )


def _format_name(name: str, times: Time, path: list[Path]) -> str:
    year, day, hour, minute = times[0, 0].to_value("yday", "date_hm").split(":")
    return name.format(year=year, day=day, hour=hour, minute=minute, stem=path[0].stem)


def _build_gsdata(
    data: np.ndarray,
    times: Time,
    freqs: un.Quantity,
    auxiliary_measurements: dict[str, np.ndarray],
    path: list[Path],
    telescope: Telescope,
    name: str,
    lst_setter: callable | None,
    history_params: dict | None = None,
    **kwargs,
) -> GSData:
    # TODO: use proper integration time...
    time_ranges = Time(
        np.hstack((times.jd, (times + 13 * un.s).jd)), format="jd"
    ).reshape((*times.shape, 2))

    if lst_setter is not None:
        kwargs["lsts"] = lst_setter(times, telescope.location)
        kwargs["lst_ranges"] = lst_setter(time_ranges, telescope.location)

    return GSData(
        data=data,
        times=times,
        freqs=freqs,
        data_unit="power",
        loads=LOADS,
        auxiliary_measurements=dict(auxiliary_measurements),
        filename=path[0] if len(path) == 1 else None,
        telescope=telescope,
        name=name,
        history=History(
            (
                Stamp(
                    "Read from ACQ file",
                    function="read_acq_to_gsdata",
                    parameters={"path": path, **(history_params or {})},
                ),
            ),
        ),
        **kwargs,
    )


def _read_selected(
    path: list[Path],
    telescope: Telescope,
    name: str,
    lst_setter: callable | None,
    selectors: dict[str, dict],
    **kwargs,
) -> GSData:
    """Read only the selected parts of a set of ACQ files.

    First the ancillary data (including times) of every file is read without
    decoding any spectra. Masks are computed from these, and then only the selected
    spectra/channels are decoded.
    """
    unknown = set(selectors) - set(SELECTORS)
    if unknown:
        raise ValueError(
            f"Unrecognized selectors: {unknown}. Available selectors: "
            f"{', '.join(SELECTORS)}"
        )

    indexed = [(pth, *_index_file(pth)) for pth in path]
    indexed = [(pth, anc, offsets) for pth, anc, offsets in indexed if len(offsets)]
    if not indexed:
        raise ACQError(f"No data in any files: {path}")

    anc = indexed[0][1]
    freqs = anc.frequencies * un.MHz
    times = time_concat(
        [Time(a.data.pop("times"), format="yday", scale="utc") for _, a, _ in indexed]
    )
    aux = {k: np.concatenate([a.data[k] for _, a, _ in indexed]) for k in anc.data}

    # Apply the selections in the same order as GSData.from_file does for readers
    # that don't select on read, so that e.g. an lst_selector indx refers to the
    # already time-selected data.
    fmask = freq_selector(freqs, **selectors.get("freq_selector", {}))
    chans = np.flatnonzero(fmask)
    if len(chans) == 0:
        raise ACQError("The freq_selector matched no frequency channels.")

    keep = np.arange(len(times))
    mask = time_selector(times, LOADS, **selectors.get("time_selector", {}))
    if mask is not None:
        keep = keep[mask]

    if "lst_selector" in selectors:
        sel_times = times[keep]
        lsts = (
            lst_setter(sel_times, telescope.location)
            if lst_setter is not None
            else sel_times.sidereal_time("apparent", telescope.location)
        )
        mask = lst_selector(lsts, LOADS, **selectors["lst_selector"])
        if mask is not None:
            keep = keep[mask]

    if len(keep) == 0:
        raise ACQError(f"Selection matched no data in files: {path}")

    # Decode the contiguous channel span, then drop any unselected channels in it.
    span = slice(chans[0], chans[-1] + 1)
    bounds = np.cumsum([0] + [len(offsets) for _, _, offsets in indexed])
    spectra = []
    for (pth, _, offsets), lo, hi in zip(indexed, bounds[:-1], bounds[1:], strict=True):
        file_keep = keep[(keep >= lo) & (keep < hi)] - lo
        if len(file_keep):
            spectra.append(
                _read_spectra(pth, offsets[file_keep], anc.meta["nfreq"], span)
            )

    data = np.concatenate(spectra, axis=1)
    if len(chans) != span.stop - span.start:
        data = data[..., chans - span.start]

    gsd = _build_gsdata(
        data=data[:, np.newaxis],
        times=times[keep],
        freqs=freqs[fmask],
        auxiliary_measurements={k: v[keep] for k, v in aux.items()},
        path=path,
        telescope=telescope,
        # Name the object by the start of the file(s), not of the selection.
        name=_format_name(name, times, path),
        lst_setter=lst_setter,
        history_params={"selectors": selectors},
        **kwargs,
    )

    if "load_selector" in selectors:
        gsd = select_loads(gsd, **selectors["load_selector"])

    return gsd


def write_gsdata_to_acq(
    gsdata: GSData,
    outfile: str | Path,
    temperature: float = 25.0,
    nblk: int = 2974,
):
    """Write a GSData object to an ACQ file.

    Parameters
    ----------
    gsdata : GSData
        The data to write.
    outfile : str or Path
        The file to write to.
    temperature : float
        The temperature of the system, in Celsius.
    nblk : int
        The number of blocks in the data.
    """
    if gsdata.nloads != 3:
        raise ValueError("Can only encode 3-load data to ACQ file.")

    ancillary = {
        "times": gsdata.times.strftime("%Y:%j:%H:%M:%S"),
        "adcmax": gsdata.auxiliary_measurements["adcmax"],
        "adcmin": gsdata.auxiliary_measurements["adcmin"],
    }

    meta = {
        "temperature": temperature,
        "nblk": nblk,
        "nfreq": gsdata.nfreqs,
        "freq_min": gsdata.freqs.min().to_value("MHz"),
        "freq_max": gsdata.freqs.max().to_value("MHz"),
        "freq_res": (gsdata.freqs[1] - gsdata.freqs[0]).to_value("MHz"),
    }

    encode(outfile, p=gsdata.data[:, 0], meta=meta, ancillary=ancillary)
