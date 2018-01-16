# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy


def xyz_to_srgb1(xyz):
    # https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
    M = numpy.array([
        [+3.2406, -1.5372, -0.4986],
        [-0.9689, +1.8758, +0.0415],
        [+0.0557, -0.2040, +1.0570],
        ])
    srgb_linear = numpy.dot(M, xyz)
    a = 0.055
    srgb = numpy.array([
        12.92 * c if c <= 0.0031308 else
        (1+a) * c**(1/2.4) - a
        for c in srgb_linear
        ])
    return srgb


def srgb1_to_xyz(srgb1):
    # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
    a = 0.055
    srgb_linear = numpy.array([
        c / 12.92 if c <= 0.04045 else
        ((c+a) / (1+a)) ** 2.4
        for c in srgb1
        ])
    M = numpy.array([
        [0.4124, 0.3576, 0.1805],
        [0.2126, 0.7152, 0.0722],
        [0.0193, 0.1192, 0.9505],
        ])
    return numpy.dot(M, srgb_linear)


def srgb1_to_srgb256(srgb1):
    return (256 * srgb1).astype(int)


def srgb256_to_srgb1(srgb256):
    return srgb256 / 256
