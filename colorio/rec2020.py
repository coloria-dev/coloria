from __future__ import division

import numpy

from .linalg import dot, solve
from .xyy import XYY


class Rec2020(object):
    """
    RGB color space from the Rec. 2020. Same as Rec. 2100, used in HDR. The
    primary colors are monochromatic at wave lengths 467nm, 532nm, and 630nm.

    <https://www.itu.int/dms_pubrec/itu-r/rec/bt/R-REC-BT.2020-2-201510-I!!PDF-E.pdf>.
    """

    def __init__(self):
        # The Y-coordinate the guessed, it does not explicitly appear in the spec.
        primaries_xyy = numpy.array(
            [[0.708, 0.292, 1 / 3], [0.170, 0.797, 1 / 3], [0.131, 0.046, 1 / 3]]
        )
        self.invM = XYY().to_xyz100(primaries_xyy.T) / 100

        self.alpha = 1.09929682680944
        self.beta = 0.018053968510807
        return

    def from_xyz100(self, xyz100):
        # TODO NaN the values smaller than 0 and larger than 1
        return solve(self.invM, xyz100 / 100)

    def to_xyz100(self, rec2020_linear):
        return 100 * dot(self.invM, rec2020_linear)

    def from_gamma(self, gamma_corrected):
        out = numpy.array(gamma_corrected, dtype=float)

        is_smaller = out <= 4.5 * self.beta
        out[is_smaller] /= 4.5
        out[~is_smaller] = ((out[~is_smaller] + self.alpha - 1) / self.alpha) ** (
            1 / 0.45
        )
        return out

    def to_gamma(self, linear):
        out = numpy.array(linear, dtype=float)

        is_smaller = linear <= self.beta
        out[is_smaller] *= 4.5
        out[~is_smaller] = self.alpha * out[~is_smaller] ** 0.45 - self.alpha + 1
        return out
