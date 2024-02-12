"""Wrapper for C-code that does the encoding-decoding."""

import ctypes
from pathlib import Path

import numpy as np

cdll = next(Path(__file__).parent.glob("decode.*.so"))
cdll = ctypes.CDLL(cdll)

_c_decode = cdll.decode
_c_decode.restype = ctypes.c_int
_c_decode.argtypes = [
    ctypes.c_char_p,
    np.ctypeslib.ndpointer(np.float64),
]

_c_encode = cdll.encode
_c_encode.restype = ctypes.c_int
_c_encode.argtypes = [
    ctypes.c_int,
    np.ctypeslib.ndpointer(np.float64),
    ctypes.c_char_p,
]


def _decode_line(line: str) -> np.ndarray:
    """
    Decode a pre-parsed line of an ACQ file.

    This is a very thin wrapper around the C `decode` function which does the same.

    Parameters
    ----------
    line : str
        A string that is uencoded. There should be no spaces at the start or end of the
        line. This is *not* verified in this function.

    Returns
    -------
    np.ndarray :
        A 1D array of float values given by the de-encoding and unpacking. These are
        arbitrarily scaled linear powers.
    """
    out = np.zeros(len(line) // 4)
    res = _c_decode(ctypes.c_char_p(line.encode("ascii")), np.ascontiguousarray(out))

    if res:
        raise SystemError("C decoder exited with an error!")

    return out


def _encode(data: np.ndarray) -> str:
    """Encode an array of data.

    This is a low-level function that takes an array of linear-scaled float data and
    converts it to a string of 64-bit encoded integers. It does *not* perform scaling
    to match the output of fastspec, but its output should be the inverse of
    :func:`_decode_line`.
    """
    out = ctypes.create_string_buffer(len(data) * 4)
    res = _c_encode(len(data), data, out)

    if res:
        raise SystemError("C encoder exited with an error!")

    return out.raw.decode("ascii")


def _encode_line(data: np.ndarray, nblk: int) -> str:
    """Encode an array of data.

    This takes an array of linear-scaled float data and converts it to a string of
    64-bit encoded integers. It performs scaling to match the output of fastspec, and
    also blanks out the first 10 entries of data (which are never used in fastspec).

    The output is *not* exactly the inverse of :func:`_decode_line`, but is the same up
    to a scaling constant, as well as some clipping on dynamic range. Thus, while it is
    important that this function scales its input (to get the data into the clipping
    range), it is *not* important to scale the output of `_decode_line`, as it is
    fully linear, and only the ratio is ever used.
    """
    # We scale the data to achieve same dynamic range as in fastspec
    d = data / (2 * nblk * len(data) * 10**3.84)

    # Set first 10 values to 0 because they are never used
    d[:10] = 0

    return _encode(d)
