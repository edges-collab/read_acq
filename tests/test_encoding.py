import pytest

import numpy as np

from read_acq import convert_file, convert_h5, decode_file, encode
from read_acq.codec import _decode_line, _encode, _encode_line


def test_roundtrip():
    a = np.logspace(-5, -1, 20)

    b = _decode_line(_encode(a))
    assert np.allclose(a, b)


def test_roundtrip_ratio():
    a = np.logspace(-5, -1, 20)
    b = np.logspace(-6, -2, 20)

    aa = _decode_line(_encode_line(a, 1000))
    bb = _decode_line(_encode_line(b, 1000))

    assert np.allclose(a[10:] / b[10:], aa[10:] / bb[10:])


@pytest.fixture(scope="module")
def sample_acq_file(tmp_path_factory):
    data = np.logspace(-5, -1, 1500).reshape((3, 5, 100))
    meta = {
        "temp": 0,
        "nblk": 1000,
        "nfreq": 100,
        "freq_min": 0,
        "freq_max": 200,
        "freq_res": 2,
    }

    fname = tmp_path_factory.mktemp("direc") / "tempfile.acq"

    ancillary = {
        "adcmin": np.ones((10, 3)) * 0.2,
        "adcmax": np.ones((10, 3)) * 0.3,
        "data_drops": np.zeros((10, 3), dtype=int),
        "times": np.array(["2016:080:01:01:01"] * 10),
    }

    encode(fname, data, meta, ancillary)
    return fname


@pytest.fixture(scope="module")
def incomplete_file(sample_acq_file):
    with open(sample_acq_file) as fl:
        lines = fl.readlines()
    last_line = lines[-1]
    meta, spec = last_line.split(" spectrum ")
    spec = _decode_line(spec.lstrip())
    spec = spec[: spec.size // 2]
    spec = _encode_line(spec, nblk=1000)
    lines[-1] = meta + " spectrum " + spec + "\n"

    new_file = sample_acq_file.parent / "incomplete.acq"
    with open(new_file, "w") as fl:
        fl.writelines(lines)
    return new_file


def test_full_file_encode(sample_acq_file):
    with open(sample_acq_file) as fl:
        assert fl.readline() == ";--temp: 0\n"


def test_full_file_h5(sample_acq_file, tmp_path_factory):
    h5_file = tmp_path_factory.mktemp("direc") / "tempfile.h5"
    convert_file(sample_acq_file, outfile=h5_file)

    new_acq = tmp_path_factory.mktemp("direc") / "tempfile_roundtrip.acq"
    # Now convert back to .acq
    convert_h5(h5_file, new_acq)

    with open(new_acq) as fl:
        assert ";--temp: 0\n" in fl.readlines()[:20]


def test_incomplete_file(incomplete_file):
    Q, p = decode_file(incomplete_file, progress=False)

    # We lost an integration!
    for pp in p:
        assert pp.T.shape == (4, 100)
    assert Q.shape == (4, 100)
