from read_acq.codec import _encode_line, _decode_line, _encode
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
