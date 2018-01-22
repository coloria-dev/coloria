# -*- coding: utf-8 -*-
#
import numpy

from .illuminants import white_point, d65
from . import cielab


class CIELCH(object):
    def __init__(self, whitepoint=white_point(d65())):
        self.cielab = cielab.CIELAB(whitepoint=whitepoint)
        return

    def from_xyz(self, xyz):
        L, u, v = self.cielab.from_xyz(xyz)
        C = numpy.sqrt(u**2 + v**2)
        h = numpy.mod(numpy.arctan2(v, u), 2*numpy.pi) / numpy.pi * 180
        return numpy.array([L, C, h])

    def to_xyz(self, lch):
        L, C, h = lch
        h_ = h * numpy.pi / 180
        lab = numpy.array([L, C * numpy.cos(h_), C * numpy.sin(h_)])
        return self.cielab.to_xyz(lab)
