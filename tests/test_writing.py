import pytest

import h5py
import numpy as np
from pathlib import Path
from scipy.io import loadmat

from read_acq import convert_file


@pytest.fixture(scope="module")
def sample_acq() -> Path:
    return Path(__file__).parent / "data/sample.acq"


@pytest.fixture(scope="module")
def tmpdir(tmp_path_factory):
    return tmp_path_factory.mktemp("direc")


def test_h5(tmpdir, sample_acq):
    outfile = tmpdir / "tempfile.h5"
    Q, p, meta = convert_file(sample_acq, outfile=outfile, write_format="h5", meta=True)

    with h5py.File(outfile, "r") as fl:
        qq = fl["spectra"]["Q"][...]

    assert np.allclose(qq[~np.isnan(qq)], Q[~np.isnan(Q)])


def test_mat(tmpdir, sample_acq):
    outfile = tmpdir / "tempfile.mat"
    Q, p, meta = convert_file(
        sample_acq, outfile=outfile, write_format="mat", meta=True
    )

    matdata = loadmat(str(outfile))

    assert np.allclose(matdata["Qratio"][~np.isnan(Q)], Q[~np.isnan(Q)])


def test_npz(tmpdir, sample_acq):
    outfile = tmpdir / "tempfile"
    Q, p, meta = convert_file(
        sample_acq, outfile=outfile, write_format="npz", meta=True
    )

    matdata = np.load(outfile.with_suffix(".npz"))

    assert np.allclose(matdata["Qratio"][~np.isnan(Q)], Q[~np.isnan(Q)])


def test_bad_writer(tmpdir, sample_acq):
    outfile = tmpdir / "tempfile"
    with pytest.raises(ValueError):
        convert_file(sample_acq, outfile=outfile, write_format="bad_format", meta=True)
