# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy


class SRGB1(object):
    def __init__(self):
        self.a = 0.055
        return

    def from_srgb_linear(self, srgb_linear):
        is_smaller = srgb_linear <= 0.0031308

        srgb = numpy.array(srgb_linear, dtype=float)
        srgb[is_smaller] *= 12.92
        srgb[~is_smaller] = (1+self.a) * srgb[~is_smaller]**(1/2.4) - self.a
        return srgb

    def to_srgb_linear(self, srgb1):
        srgb_linear = numpy.array(srgb1, dtype=float)

        # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
        is_smaller = srgb_linear <= 12.92 * 0.0031308  # 0.040449936

        srgb_linear[is_smaller] /= 12.92
        srgb_linear[~is_smaller] = (
                (srgb_linear[~is_smaller] + self.a) / (1+self.a)
                )**2.4
        return srgb_linear

    def srgb_gamut(self, filename='srgb.vtu', n=50):
        import meshio
        import meshzoo
        points, cells = meshzoo.cube(nx=n, ny=n, nz=n)
        meshio.write(
            filename,
            points, {'tetra': cells},
            point_data={'srgb1': self.from_srgb_linear(points)}
            )
        return
