"""Package providing functionality for reading .ACQ format files."""

__all__ = ["decode_file", "encode", "read_acq_to_gsdata", "write_gsdata_to_acq"]

from importlib.metadata import PackageNotFoundError, version

from .gsdata import read_acq_to_gsdata, write_gsdata_to_acq
from .read_acq import decode_file, encode

try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: no cover
    # package is not installed
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
