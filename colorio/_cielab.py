import numpy

from ._color_space import ColorSpace
from .illuminants import whitepoints_cie1931


class CIELAB(ColorSpace):
    def __init__(self, whitepoint=whitepoints_cie1931["D65"]):
        self.whitepoint = whitepoint
        self.labels = ["L*", "a*", "b*"]
        return

    def from_xyz100(self, xyz):
        def f(t):
            delta = 6 / 29
            out = numpy.array(t, dtype=float)
            is_greater = out > delta ** 3
            out[is_greater] = 116 * numpy.cbrt(out[is_greater]) - 16
            out[~is_greater] = out[~is_greater] / (delta / 2) ** 3
            return out

        fx, fy, fz = f((xyz.T / self.whitepoint).T)
        return numpy.array([fy, 125 / 29 * (fx - fy), 50 / 29 * (fy - fz)])

    def to_xyz100(self, lab):
        def f1(t):
            delta = 6.0 / 29.0
            out = numpy.array(t, dtype=float)
            is_greater = out > 2.0 / 29.0
            out[is_greater] = (out[is_greater] + 4.0 / 29.0) ** 3
            out[~is_greater] = 3 * delta ** 2 * out[~is_greater]
            return out

        L, a, b = lab
        return (
            f1(numpy.array([L / 116 + a / 500, L / 116, L / 116 - b / 200])).T
            * self.whitepoint
        ).T
