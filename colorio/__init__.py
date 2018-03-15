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

from .cam16 import CAM16, CAM16UCS
from .ciecam02 import CIECAM02, CAM02, NegativeAError
from .ciehcl import CIEHCL
from .cielab import CIELAB
from .cielch import CIELCH
from .cieluv import CIELUV
from . import illuminants
from .ictcp import ICtCp
from .ipt import IPT
from .jzazbz import JzAzBz
from . import observers
from .rec2020 import Rec2020
from .srgb import SrgbLinear
from .xyy import XYY
from .xyz import XYZ
from .tools import (
    show_srgb_gamut, show_hdr_gamut, show_visible_gamut, show_flat_gamut,
    delta, show_ebner_fairchild, show_hung_berns, show_munsell, show_macadam,
    plot_macadam, show_luo_rigg, show_xiao, show_straights
    )

try:
    import pipdate
except ImportError:
    pass
else:
    if pipdate.needs_checking(__name__):
        print(pipdate.check(__name__, __version__), end='')
