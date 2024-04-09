"""Package providing functionality for reading .ACQ format files."""

import contextlib

__all__ = ["decode_file", "encode", "read_acq_to_gsdata", "write_gsdata_to_acq"]
from pkg_resources import DistributionNotFound, get_distribution

from .gsdata import read_acq_to_gsdata, write_gsdata_to_acq
from .read_acq import decode_file, encode

with contextlib.suppress(DistributionNotFound):
    __version__ = get_distribution(__name__).version

del get_distribution
del DistributionNotFound
del contextlib
