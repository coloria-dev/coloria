# -*- coding: utf-8 -*-
#
import numpy

from .illuminants import white_point, d65
from . import cielab


def from_xyz(xyz, whitepoint=white_point(d65())):
    L, u, v = cielab.from_xyz(xyz, whitepoint=whitepoint)
    C = numpy.sqrt(u**2 + v**2)
    h = numpy.arctan2(v, u)
    return numpy.array([L, C, h])


def to_xyz(lch, whitepoint=white_point(d65())):
    L, C, h = lch
    lab = numpy.array([L, C * numpy.cos(h), C * numpy.sin(h)])
    return cielab.to_xyz(lab, whitepoint=whitepoint)
