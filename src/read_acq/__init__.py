"""Package providing functionality for reading .ACQ format files."""
from pkg_resources import DistributionNotFound, get_distribution

from . import codec, writers
from .read_acq import (
    Ancillary,
    convert_file,
    convert_h5,
    decode_file,
    decode_files,
    encode,
)

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

del get_distribution
del DistributionNotFound
