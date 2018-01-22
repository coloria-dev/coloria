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

    def from_xyz100(self, xyz):
        # https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
        # http://www.color.org/srgb.pdf
        # TODO NaN the values smaller than 0 and larger than 1
        return numpy.einsum('ij,j...->i...', self.M, xyz) / 100

    def to_xyz100(self, srgb1_linear):
        # https://stackoverflow.com/a/48387507/353337
        x = numpy.array(srgb1_linear)
        return 100 * numpy.linalg.solve(
            self.M, x.reshape(x.shape[0], -1)
            ).reshape(x.shape)

    # pylint: disable=no-self-use
    def from_srgb1(self, srgb1):
        srgb_linear = numpy.array(srgb1, dtype=float)

        a = 0.055
        # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
        is_smaller = srgb_linear <= 12.92 * 0.0031308  # 0.040449936

        srgb_linear[is_smaller] /= 12.92
        srgb_linear[~is_smaller] = ((srgb_linear[~is_smaller] + a)/(1+a))**2.4
        return srgb_linear

    def to_srgb1(self, srgb_linear):
        a = 0.055
        is_smaller = srgb_linear <= 0.0031308
        srgb = numpy.array(srgb_linear, dtype=float)
        srgb[is_smaller] *= 12.92
        srgb[~is_smaller] = (1+a) * srgb[~is_smaller]**(1/2.4) - a
        return srgb
