import numpy

from ._cieluv import CIELUV
from ._color_space import ColorSpace
from .illuminants import whitepoints_cie1931


class CIEHCL(ColorSpace):
    def __init__(self, whitepoint=whitepoints_cie1931["D65"]):
        super().__init__()
        self.cieluv = CIELUV(whitepoint=whitepoint)
        self.labels = ["L", "C", "h"]
        self.k0 = 0

    def from_xyz100(self, xyz):
        L, u, v = self.cieluv.from_xyz100(xyz)
        C = numpy.hypot(u, v)
        h = numpy.mod(numpy.arctan2(v, u), 2 * numpy.pi) / numpy.pi * 180
        return numpy.array([L, C, h])

    def to_xyz100(self, lch):
        L, C, h = lch
        h_ = h * numpy.pi / 180
        luv = numpy.array([L, C * numpy.cos(h_), C * numpy.sin(h_)])
        return self.cieluv.to_xyz100(luv)
