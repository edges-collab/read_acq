"""Tests of writing the data."""
from pathlib import Path

import h5py
import numpy as np
import pytest
from read_acq import convert_file


@pytest.fixture(scope="module")
def sample_acq() -> Path:
    return Path(__file__).parent / "data/sample.acq"


@pytest.fixture(scope="module")
def tmpdir(tmp_path_factory):
    return tmp_path_factory.mktemp("direc")


def test_h5(tmpdir, sample_acq):
    outfile = tmpdir / "tempfile.h5"
    q, p, meta = convert_file(sample_acq, outfile=outfile, write_format="h5", meta=True)

    with h5py.File(outfile, "r") as fl:
        qq = fl["spectra"]["Q"][...]

    assert np.allclose(qq[~np.isnan(qq)], q[~np.isnan(q)])


def test_bad_writer(tmpdir, sample_acq):
    outfile = tmpdir / "tempfile"
    with pytest.raises(ValueError, match="does not have an associated writer"):
        convert_file(sample_acq, outfile=outfile, write_format="bad_format", meta=True)
