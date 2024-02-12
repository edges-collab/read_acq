"""Package providing functionality for reading .ACQ format files."""

import contextlib

__all__ = ["convert_file", "decode_file", "encode"]
from pkg_resources import DistributionNotFound, get_distribution

from .read_acq import convert_file, decode_file, encode

with contextlib.suppress(DistributionNotFound):
    __version__ = get_distribution(__name__).version

del get_distribution
del DistributionNotFound
del contextlib
