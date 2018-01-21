# -*- coding: utf-8 -*-
#
from __future__ import print_function

from .__about__ import (
    __author__,
    __email__,
    __copyright__,
    __credits__,
    __license__,
    __version__,
    __maintainer__,
    __status__,
    )

from .cam02 import CAM02_UCS
from .ciecam02 import CIECAM02
from .ciehcl import CIEHCL
from .cielab import CIELAB
from .cielch import CIELCH
from .cieluv import CIELUV
from . import illuminants
from .jzazbz import JzAzBz
from . import observers
from .srgb import SRGB1, SrgbLinear
from .xyy import XYY
from . import xyz

# try:
#     import pipdate
# except ImportError:
#     pass
# else:
#     if pipdate.needs_checking(__name__):
#         print(pipdate.check(__name__, __version__))
