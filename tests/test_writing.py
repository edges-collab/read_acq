from pathlib import Path
from read_acq import decode_file
import h5py
import numpy as np
from scipy.io import loadmat


def test_h5(tmp_path_factory):
    data = Path(__file__).parent / "data/sample.acq"

    outfile = tmp_path_factory.mktemp("direc") / "tempfile.h5"
    Q, p, meta = decode_file(data, outfile=outfile, write_formats=["h5"], meta=True)

    with h5py.File(outfile, "r") as fl:
        qq = fl["spectra"]["Q"][...]

    assert np.allclose(qq[~np.isnan(qq)], Q[~np.isnan(Q)])


def test_mat(tmp_path_factory):
    data = Path(__file__).parent / "data/sample.acq"
    outfile = tmp_path_factory.mktemp("direc") / "tempfile.mat"
    Q, p, meta = decode_file(data, outfile=outfile, write_formats=["mat"], meta=True)

    matdata = loadmat(outfile)

    assert np.allclose(matdata["Qratio"][~np.isnan(Q)], Q[~np.isnan(Q)])


def test_npz(tmp_path_factory):
    data = Path(__file__).parent / "data/sample.acq"
    outfile = tmp_path_factory.mktemp("direc") / "tempfile.npz"
    Q, p, meta = decode_file(data, outfile=outfile, write_formats=["npz"], meta=True)

    matdata = np.load(outfile)

    assert np.allclose(matdata["Qratio"][~np.isnan(Q)], Q[~np.isnan(Q)])
