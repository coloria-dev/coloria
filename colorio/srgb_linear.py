# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy


class SrgbLinear(object):
    def __init__(self):
        self.M = numpy.array([
            [+3.2406255, -1.537208, -0.4986286],
            [-0.9689307, +1.8757561, +0.0415175],
            [+0.0557101, -0.2040211, +1.0569959],
            ])
        return

    def from_xyz(self, xyz):
        # https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
        # http://www.color.org/srgb.pdf
        return numpy.dot(self.M, xyz)

    def to_xyz(self, srgb1_linear):
        return numpy.linalg.solve(self.M, srgb1_linear)
