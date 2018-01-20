# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .illuminants import white_point, d65
from . import srgb_linear
from . import srgb1


def from_xyz(xyz, whitepoint=white_point(d65())):
    def f(t):
        delta = 6.0/29.0
        out = numpy.array(t, dtype=float)
        is_greater = out > delta**3
        out[is_greater] = numpy.cbrt(out[is_greater])
        out[~is_greater] = out[~is_greater]/3/delta**2 + 4.0/29.0
        return out

    fx, fy, fz = f((xyz.T / whitepoint).T)
    return numpy.array([
        116 * fy - 16,
        500 * (fx - fy),
        200 * (fy - fz),
        ])


def to_xyz(cielab, whitepoint=white_point(d65())):
    def f1(t):
        delta = 6.0/29.0
        out = numpy.array(t, dtype=float)
        is_greater = out > delta
        out[is_greater] = out[is_greater]**3
        out[~is_greater] = 3*delta**2 * (out[~is_greater] - 4.0/29.0)
        return out

    L, a, b = cielab
    return (f1(numpy.array([
        (L+16)/116 + a/500,
        (L+16)/116,
        (L+16)/116 - b/200,
        ])).T * whitepoint).T


def srgb_gamut(filename='srgb-cielab.vtu', n=50):
    import meshio
    import meshzoo
    points, cells = meshzoo.cube(nx=n, ny=n, nz=n)
    pts = from_xyz(srgb_linear.to_xyz(points.T)).T
    rgb = srgb1.from_srgb_linear(points)
    meshio.write(
        filename,
        pts, {'tetra': cells},
        point_data={'srgb': rgb}
        )
    return
