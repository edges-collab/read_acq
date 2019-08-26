from os import path
import numpy as np
from scipy import io
import copy

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
        print("Writing {}...".format(outfile), end='', flush=True)
        func(outfile, ancillary, **data)
        print(" done")

    return writer_wrapper


def _get_fname(outfile=None, format=None):
    if not outfile:
        raise Exception("You need to provide either outfile or infile!")

    if path.splitext(outfile)[1] in ['.mat', '.h5', '.npz', '.acq']:
        outfile = path.splitext(outfile)[0] + "." + format
    else:
        outfile = outfile + "." + format

    return outfile


@writer
def _write_h5(outfile=None, ancillary=None, **data):
    if not HAVE_H5PY:
        raise RuntimeError("You need the h5py library to write to HDF5 format!")

    with h5py.File(outfile, 'w') as fl:
        if ancillary:
            for k, v in ancillary.items():
                fl.attrs[k] = v

        for k, v in data.items():
            fl[k] = v


@writer
def _write_npz(outfile=None, ancillary=None, **data):
    np.savez(outfile, **data, **ancillary)


@writer
def _write_mat(outfile=None, ancillary=None, **data):
    data = copy.deepcopy(data)
    data.update(ancillary)
    io.savemat(outfile, data)