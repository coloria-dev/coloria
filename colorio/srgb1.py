# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy


def from_srgb_linear(srgb_linear):
    a = 0.055
    is_smaller = srgb_linear <= 0.0031308

    srgb = numpy.array(srgb_linear, dtype=float)
    srgb[is_smaller] *= 12.92
    srgb[~is_smaller] = (1+a) * srgb[~is_smaller]**(1/2.4) - a
    return srgb


def to_srgb_linear(srgb1):
    srgb_linear = numpy.array(srgb1, dtype=float)

    # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
    is_smaller = srgb_linear <= 12.92 * 0.0031308  # 0.040449936

    a = 0.055
    srgb_linear[is_smaller] /= 12.92
    srgb_linear[~is_smaller] = ((srgb_linear[~is_smaller] + a) / (1+a))**2.4
    return srgb_linear


def to_srgb256(srgb1):
    return (256 * srgb1).astype(int)


def from_srgb256(srgb256):
    return srgb256 / 256


def srgb_gamut(filename='srgb.vtu', n=50):
    import meshio
    import meshzoo
    points, cells = meshzoo.cube(nx=n, ny=n, nz=n)
    meshio.write(
        filename,
        points, {'tetra': cells},
        point_data={'srgb1': from_srgb_linear(points)}
        )
    return
