# -*- coding: utf-8 -*-
#
import numpy

from .illuminants import white_point, d65
from . import cieluv


class CIEHCL(object):

    def __init__(self, whitepoint=white_point(d65())):
        self.cieluv = cieluv.CIELUV(whitepoint=whitepoint)
        return

    def from_xyz(self, xyz):
        L, u, v = self.cieluv.from_xyz(xyz)
        C = numpy.sqrt(u**2 + v**2)
        h = numpy.arctan2(v, u)
        return numpy.array([L, C, h])

    def to_xyz(self, lch):
        L, C, h = lch
        luv = numpy.array([L, C * numpy.cos(h), C * numpy.sin(h)])
        return self.cieluv.to_xyz(luv)
