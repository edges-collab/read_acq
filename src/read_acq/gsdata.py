"""Read an ACQ file into a GSData object."""
from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
from astropy import units as un
from astropy.coordinates import EarthLocation, UnknownSiteException
from astropy.time import Time
from pygsdata import GSData
from pygsdata.constants import KNOWN_LOCATIONS

from .read_acq import ACQError, decode_file


def read_acq_to_gsdata(
    path: str | Path | Sequence[str | Path],
    telescope_location: str | EarthLocation,
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

    if isinstance(telescope_location, str):
        try:
            telescope_location = EarthLocation.of_site(telescope_location)
        except UnknownSiteException as e:
            try:
                telescope_location = KNOWN_LOCATIONS[telescope_location]
            except KeyError:
                raise ValueError(
                    "telescope_location must be an EarthLocation or a known site, "
                    f"got {telescope_location}. Known sites are {KNOWN_LOCATIONS}."
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
        times = np.concatenate((times, times_))
        anc.data = {k: np.concatenate((v, anc_.data[k])) for k, v in anc.data.items()}

    if pant.size == 0:
        raise ACQError(f"No data in any files: {path}")

    year, day, hour, minute = times[0, 0].to_value("yday", "date_hm").split(":")
    name = name.format(year=year, day=day, hour=hour, minute=minute, stem=path[0].stem)
    return GSData(
        data=np.array([pant.T, pload.T, plns.T])[:, np.newaxis],
        time_array=times,
        freq_array=anc.frequencies * un.MHz,
        data_unit="power",
        loads=("ant", "internal_load", "internal_load_plus_noise_source"),
        auxiliary_measurements={name: anc.data[name] for name in anc.data},
        filename=path[0] if len(path) == 1 else None,
        telescope_location=telescope_location,
        name=name,
        **kwargs,
    )
