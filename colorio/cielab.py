# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy


def xyz_to_cielab(xyz, whitepoint):
    def f(t):
        delta = 6.0/29.0
        out = numpy.array(t, dtype=float)
        is_greater = out > delta**3
        is_smaller = numpy.logical_not(is_greater)
        out[is_greater] = numpy.cbrt(out[is_greater])
        out[is_smaller] = out[is_smaller]/3/delta**2 + 4.0/29.0
        return out

    fxyz = f(xyz / whitepoint)
    return numpy.array([
        116 * fxyz[1] - 16,
        500 * (fxyz[0] - fxyz[1]),
        200 * (fxyz[1] - fxyz[2]),
        ])


def cielab_to_xyz(cielab, whitepoint):
    def f1(t):
        delta = 6.0/29.0
        out = numpy.array(t, dtype=float)
        is_greater = out > delta
        is_smaller = numpy.logical_not(is_greater)
        out[is_greater] = out[is_greater]**3
        out[is_smaller] = 3*delta**2 + (out[is_smaller] - 4.0/29.0)

    return whitepoint * f1(numpy.array([
        (cielab[0]+16)/116 + cielab[1]/500,
        (cielab[0]+16)/116,
        (cielab[0]+16)/116 - cielab[2]/500,
        ]))
