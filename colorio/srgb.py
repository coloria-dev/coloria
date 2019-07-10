from __future__ import division

import numpy

from .illuminants import whitepoints_cie1931
from .xyy import XYY
from .linalg import dot, solve


class SrgbLinear(object):
    def __init__(self, whitepoint_correction=True):
        # The standard actually gives the values in terms of M, but really inv(M) is a
        # direct derivative of the primary specification at
        # <https://en.wikipedia.org/wiki/SRGB>.
        primaries_xyy = numpy.array(
            [[0.64, 0.33, 0.2126], [0.30, 0.60, 0.7152], [0.15, 0.06, 0.0722]]
        )
        self.invM = XYY().to_xyz100(primaries_xyy.T)

        if whitepoint_correction:
            # The above values are given only approximately, resulting in the fact that
            # SRGB(1.0, 1.0, 1.0) is only approximately mapped into the reference
            # whitepoint D65. Add a correction here.
            correction = whitepoints_cie1931["D65"] / numpy.sum(self.invM, axis=1)
            self.invM = (self.invM.T * correction).T

        self.invM /= 100
        # numpy.linalg.inv(self.invM) is the matrix in the spec:
        # M = numpy.array([
        #     [+3.2406255, -1.537208, -0.4986286],
        #     [-0.9689307, +1.8757561, +0.0415175],
        #     [+0.0557101, -0.2040211, +1.0569959],
        # ])
        # self.invM = numpy.linalg.inv(M)
        self.labels = ["R", "G", "B"]
        return

    def from_xyz100(self, xyz):
        # https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
        # http://www.color.org/srgb.pdf
        # TODO NaN the values smaller than 0 and larger than 1
        return solve(self.invM, xyz / 100)

    def to_xyz100(self, srgb1_linear):
        # Note: The Y value is often used for grayscale conversion.
        # 0.2126 * R_linear + 0.7152 * G_linear + 0.0722 * B_linear
        return 100 * dot(self.invM, srgb1_linear)

    def from_srgb1(self, srgb1):
        srgb_linear = numpy.array(srgb1, dtype=float)

        a = 0.055
        # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
        is_smaller = srgb_linear <= 0.040449936  # 12.92 * 0.0031308

        srgb_linear[is_smaller] /= 12.92
        srgb_linear[~is_smaller] = ((srgb_linear[~is_smaller] + a) / (1 + a)) ** 2.4
        return srgb_linear

    def to_srgb1(self, srgb_linear):
        a = 0.055
        is_smaller = srgb_linear <= 0.0031308
        srgb = numpy.array(srgb_linear, dtype=float)
        srgb[is_smaller] *= 12.92
        srgb[~is_smaller] = (1 + a) * srgb[~is_smaller] ** (1 / 2.4) - a
        return srgb
