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
    return numpy.dot(_M, xyz)


def to_xyz(srgb1_linear):
    return numpy.linalg.solve(_M, srgb1_linear)
