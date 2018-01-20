# -*- coding: utf-8 -*-
#
import numpy

from .illuminants import white_point, d65
from . import cieluv


def from_xyz(xyz, whitepoint=white_point(d65())):
    L, u, v = cieluv.from_xyz(xyz, whitepoint=whitepoint)
    C = numpy.sqrt(u**2 + v**2)
    h = numpy.arctan2(v, u)
    return numpy.array([L, C, h])


def to_xyz(lch, whitepoint=white_point(d65())):
    L, C, h = lch
    luv = numpy.array([L, C * numpy.cos(h), C * numpy.sin(h)])
    return cieluv.to_xyz(luv, whitepoint=whitepoint)
