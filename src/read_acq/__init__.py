from pkg_resources import get_distribution, DistributionNotFound
from .read_acq import (
    Ancillary,
    decode_files,
    decode_file,
    encode,
    convert_h5,
    convert_file,
)
from . import writers
from . import codec

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

del get_distribution
del DistributionNotFound
