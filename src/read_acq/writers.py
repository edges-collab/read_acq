from os import path
import numpy as np
from scipy import io
import copy
from edges_io.h5 import HDF5RawSpectrum

try:
    import h5py

    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

_WRITERS = []


def writer(func):
    format = func.__name__.split("_")[-1]
    _WRITERS.append(format)

    def writer_wrapper(outfile=None, ancillary=None, **data):
        outfile = _get_fname(outfile, format)
        print(f"Writing {outfile}...", end="", flush=True)
        func(outfile, ancillary, **data)
        print(" done")

    return writer_wrapper


def _get_fname(outfile=None, format=None):
    if not outfile:
        raise Exception("You need to provide either outfile or infile!")

    if path.splitext(outfile)[1] in [".mat", ".h5", ".npz", ".acq"]:
        outfile = path.splitext(outfile)[0] + "." + format
    else:
        outfile = outfile + "." + format

    return outfile


@writer
def _write_h5(outfile=None, ancillary=None, **data):
    """Writes HDF5 file exactly as FastSpec should now write it out."""
    if not HAVE_H5PY:
        raise RuntimeError("You need the h5py library to write to HDF5 format!")

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
