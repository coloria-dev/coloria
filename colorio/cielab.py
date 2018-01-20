# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .illuminants import white_point, d65
from .srgb_linear import SrgbLinear
from .srgb1 import SRGB1


class CIELAB(object):
    def __init__(self, whitepoint=white_point(d65())):
        self.whitepoint = whitepoint
        return

    def from_xyz(self, xyz):
        def f(t):
            delta = 6.0/29.0
            out = numpy.array(t, dtype=float)
            is_greater = out > delta**3
            out[is_greater] = numpy.cbrt(out[is_greater])
            out[~is_greater] = out[~is_greater]/3/delta**2 + 4.0/29.0
            return out

        fx, fy, fz = f((xyz.T / self.whitepoint).T)
        return numpy.array([
            116 * fy - 16,
            500 * (fx - fy),
            200 * (fy - fz),
            ])

    def to_xyz(self, lab):
        def f1(t):
            delta = 6.0/29.0
            out = numpy.array(t, dtype=float)
            is_greater = out > delta
            out[is_greater] = out[is_greater]**3
            out[~is_greater] = 3*delta**2 * (out[~is_greater] - 4.0/29.0)
            return out

        L, a, b = lab
        return (f1(numpy.array([
            (L+16)/116 + a/500,
            (L+16)/116,
            (L+16)/116 - b/200,
            ])).T * self.whitepoint).T

    def srgb_gamut(self, filename='srgb-cielab.vtu', n=50):
        import meshio
        import meshzoo
        points, cells = meshzoo.cube(nx=n, ny=n, nz=n)
        pts = self.from_xyz(SrgbLinear().to_xyz(points.T)).T
        rgb = SRGB1().from_srgb_linear(points)
        meshio.write(
            filename,
            pts, {'tetra': cells},
            point_data={'srgb': rgb}
            )
        return
