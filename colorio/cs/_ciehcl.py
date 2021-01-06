import numpy

from ..illuminants import whitepoints_cie1931
from ._cieluv import CIELUV
from ._color_space import ColorSpace


class CIEHCL(ColorSpace):
    def __init__(self, whitepoint=whitepoints_cie1931["D65"]):
        super().__init__("CIEHCL", ["L", "C", "h"], 0)
        self.cieluv = CIELUV(whitepoint=whitepoint)

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
