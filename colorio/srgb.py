# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .xyy import XYY
from .linalg import dot, solve


class SrgbLinear(object):
    def __init__(self):
        # The standard actually gives the values in terms of M, but really
        # inv(M) is a direct derivative of the primary specification.
        primaries_xyy = numpy.array(
            [[0.64, 0.33, 0.2126], [0.30, 0.60, 0.7152], [0.15, 0.06, 0.0722]]
        )
        self.invM = XYY().to_xyz100(primaries_xyy.T) / 100
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
        return 100 * dot(self.invM, srgb1_linear)

    # pylint: disable=no-self-use
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
