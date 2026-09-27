"""Wrapper for C-code that does the encoding-decoding."""

import ctypes
from importlib.resources import as_file, files

import numpy as np

# Find the compiled library through importlib.resources rather than __file__.
_lib = min(
    (p for p in files(__package__).iterdir() if p.name.startswith("libdecode.")),
    key=lambda p: p.name,
)
with as_file(_lib) as _path:
    cdll = ctypes.CDLL(str(_path))

_c_decode = cdll.decode
_c_decode.restype = ctypes.c_int
_c_decode.argtypes = [
    ctypes.c_char_p,
    np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
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
    _decode_into(line.encode("ascii"), out)
    return out


def _decode_into(line: bytes, out: np.ndarray) -> None:
    """
    Decode a pre-parsed, uencoded byte-string into a pre-allocated array.

    Each group of four characters decodes independently to one value, so any
    four-character-aligned substring of a spectrum may be passed to decode just those
    channels.

    Parameters
    ----------
    line : bytes
        The uencoded data, with no leading spaces.
    out : np.ndarray
        A contiguous float64 array of length at least ``len(line) // 4``, into which
        the decoded values are written.
    """
    if _c_decode(ctypes.c_char_p(line), out) > 0:
        raise SystemError("C decoder exited with an error!")


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
