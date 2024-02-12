"""Functions that write ACQ file data to different formats."""
from __future__ import annotations

from pathlib import Path

_WRITERS = []


def writer(func):
    """Convert a function to a writer function."""
    fmt = func.__name__.split("_")[-1]
    _WRITERS.append(fmt)

    def writer_wrapper(outfile=None, ancillary=None, **data):
        outfile = _get_fname(outfile, fmt)
        func(outfile, ancillary, **data)

    return writer_wrapper


def _get_fname(outfile: Path | str | None = None, fmt: str | None = None):
    if not outfile:
        raise Exception("You need to provide either outfile or infile!")
    outfile = Path(outfile)
    return outfile.with_suffix(f".{fmt}")


@writer
def _write_h5(outfile=None, ancillary=None, **data):
    """Write a HDF5 file exactly as FastSpec should now write it out."""
    try:
        from edges_io.h5 import HDF5RawSpectrum
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "To write to h5, you need to install edges_io or do "
            "`pip install read_acq[h5]"
        ) from e

    spectra = {
        "Q": data["Qratio"],
        "p0": data["p0"],
        "p1": data["p1"],
        "p2": data["p2"],
    }
    ancillary.update(fastspec_version=data["fastspec_version"])
    freq_anc = {"frequencies": data["freqs"]}
    time_anc = {name: data["time_data"][name] for name in data["time_data"]}

    obj = HDF5RawSpectrum.from_data(
        {
            "spectra": spectra,
            "meta": ancillary,
            "freq_ancillary": freq_anc,
            "time_ancillary": time_anc,
        }
    )

    obj.write(outfile)
