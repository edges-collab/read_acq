from read_acq.codec import _encode_line, _decode_line, _encode
from read_acq import encode

import numpy as np


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


def test_full_file_encode(tmp_path_factory):
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
        "time": np.array(["time"] * 10),
    }

    encode(fname, data, meta, ancillary)

    with open(fname) as fl:
        assert fl.readline() == ";--temp: 0\n"
