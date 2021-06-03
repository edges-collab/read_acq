"""Functions that write ACQ file data to different formats."""
import copy
from os import path

import numpy as np
from scipy import io

_WRITERS = []


def writer(func):
    """Decorator to denote a function as a writer."""
    fmt = func.__name__.split("_")[-1]
    _WRITERS.append(fmt)

    def writer_wrapper(outfile=None, ancillary=None, **data):
        outfile = _get_fname(outfile, fmt)
        print(f"Writing {outfile}...", end="", flush=True)
        func(outfile, ancillary, **data)
        print(" done")

    return writer_wrapper


def _get_fname(outfile=None, fmt=None):
    if not outfile:
        raise Exception("You need to provide either outfile or infile!")

    if path.splitext(outfile)[1] in [".mat", ".h5", ".npz", ".acq"]:
        outfile = path.splitext(outfile)[0] + "." + fmt
    else:
        outfile = outfile + "." + fmt

    return outfile


@writer
def _write_h5(outfile=None, ancillary=None, **data):
    """Writes HDF5 file exactly as FastSpec should now write it out."""
    try:
        from edges_io.h5 import HDF5RawSpectrum
    except ImportError:
        raise ImportError(
            "To write to h5, you need to install edges_io or do "
            "`pip install read_acq[h5]"
        )

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


@writer
def _write_npz(outfile=None, ancillary=None, **data):
    data = copy.deepcopy(data)
    data.update(ancillary)
    np.savez(outfile, **data)


@writer
def _write_mat(outfile=None, ancillary=None, **data):
    data = copy.deepcopy(data)
    data.update(ancillary)
    io.savemat(outfile, data)
