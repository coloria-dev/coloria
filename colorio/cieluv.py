# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .illuminants import whitepoints_cie1931


class CIELUV(object):
    def __init__(self, whitepoint=whitepoints_cie1931['D65']):
        self.whitepoint = whitepoint
        return

    def from_xyz100(self, xyz):
        def f(t):
            delta = 6.0/29.0
            out = numpy.array(t, dtype=float)
            is_greater = out > delta**3
            out[is_greater] = 116*numpy.cbrt(out[is_greater]) - 16
            out[~is_greater] = out[~is_greater] / (delta/2)**3
            return out

        L = f(xyz[1]/self.whitepoint[1])

        x, y, z = xyz
        u = 4*x / (x + 15*y + 3*z)
        v = 9*y / (x + 15*y + 3*z)
        wx, wy, wz = self.whitepoint
        un = 4*wx / (wx + 15*wy + 3*wz)
        vn = 9*wy / (wx + 15*wy + 3*wz)
        return numpy.array([L, 13*L*(u - un), 13*L*(v - vn)])

    def to_xyz100(self, luv):
        def f1(t):
            out = numpy.array(t, dtype=float)
            is_greater = out > 8
            out[is_greater] = ((out[is_greater]+16)/116)**3
            out[~is_greater] = out[~is_greater] * (3.0/29.0)**3
            return out

        L, u, v = luv

        wx, wy, wz = self.whitepoint
        un = 4*wx / (wx + 15*wy + 3*wz)
        vn = 9*wy / (wx + 15*wy + 3*wz)

        uu = u/(13*L) + un
        vv = v/(13*L) + vn

        Y = wy * f1(L)
        X = Y * 9*uu/(4*vv)
        Z = Y * (12 - 3*uu - 20*vv) / (4*vv)
        return numpy.array([X, Y, Z])
