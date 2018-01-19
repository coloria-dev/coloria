# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

_M = numpy.array([
    [+3.2406255, -1.537208, -0.4986286],
    [-0.9689307, +1.8757561, +0.0415175],
    [+0.0557101, -0.2040211, +1.0569959],
    ])


def from_xyz(xyz):
    # https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
    # http://www.color.org/srgb.pdf
    srgb_linear = numpy.dot(_M, xyz)

    a = 0.055
    is_smaller = srgb_linear <= 0.0031308
    is_greater = numpy.logical_not(is_smaller)

    srgb = srgb_linear
    srgb[is_smaller] *= 12.92
    srgb[is_greater] = (1+a) * srgb[is_greater]**(1/2.4) - a
    return srgb


def to_xyz(srgb1):
    srgb_linear = numpy.array(srgb1, dtype=float)

    # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
    is_smaller = srgb_linear <= 0.04045
    is_greater = numpy.logical_not(is_smaller)

    a = 0.055
    srgb_linear[is_smaller] /= 12.92
    srgb_linear[is_greater] = ((srgb_linear[is_greater] + a) / (1+a))**2.4

    return numpy.linalg.solve(_M, srgb_linear)


def to_srgb256(srgb1):
    return (256 * srgb1).astype(int)


def from_srgb256(srgb256):
    return srgb256 / 256
