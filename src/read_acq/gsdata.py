"""Read an ACQ file into a GSData object."""
from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
from astropy import units as un
from astropy.time import Time
from pygsdata import KNOWN_TELESCOPES, GSData, Telescope
from pygsdata.readers import gsdata_reader
from pygsdata.utils import time_concat

from .read_acq import ACQError, decode_file, encode


@gsdata_reader(select_on_read=False, formats=["acq"])
def read_acq_to_gsdata(
    path: str | Path | Sequence[str | Path],
    telescope: Telescope = KNOWN_TELESCOPES["edges-low"],
    name: str = "{year}-{day}:{hour}:{minute}",
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
    **kwargs
        Additional keyword arguments to pass to the GSData constructor.
    """
    path = [Path(path)] if isinstance(path, (str, Path)) else [Path(p) for p in path]

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

    # Read the _first_ file to get the metadata
    _, (pant, pload, plns), anc = decode_file(path[0], meta=True)
    times = Time(anc.data.pop("times"), format="yday", scale="utc")

    # Concatenate all the files
    for p in path[1:]:
        _, (pant_, pload_, plns_), anc_ = decode_file(p, meta=True)
        times_ = Time(anc_.data.pop("times"), format="yday", scale="utc")

        pant = np.concatenate((pant, pant_), axis=1)
        pload = np.concatenate((pload, pload_), axis=1)
        plns = np.concatenate((plns, plns_), axis=1)
        times = time_concat((times, times_))
        anc.data = {k: np.concatenate((v, anc_.data[k])) for k, v in anc.data.items()}

    if pant.size == 0:
        raise ACQError(f"No data in any files: {path}")

    year, day, hour, minute = times[0, 0].to_value("yday", "date_hm").split(":")
    name = name.format(year=year, day=day, hour=hour, minute=minute, stem=path[0].stem)
    return GSData(
        data=np.array([pant.T, pload.T, plns.T])[:, np.newaxis],
        times=times,
        freqs=anc.frequencies * un.MHz,
        data_unit="power",
        loads=("ant", "internal_load", "internal_load_plus_noise_source"),
        auxiliary_measurements={name: anc.data[name] for name in anc.data},
        filename=path[0] if len(path) == 1 else None,
        telescope=telescope,
        name=name,
        **kwargs,
    )


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
