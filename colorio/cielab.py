# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .illuminants import whitepoints_cie1931


class CIELAB(object):
    def __init__(self, whitepoint=whitepoints_cie1931['D65']):
        self.whitepoint = whitepoint
        self.labels = ['L*', 'a*', 'b*']
        return

    def from_xyz100(self, xyz):
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

    def to_xyz100(self, lab):
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
