# -*- coding: utf-8 -*-
#
import numpy

from .illuminants import white_point, d65
from . import cieluv


class CIEHCL(object):

    def __init__(self, whitepoint=white_point(d65())):
        self.cieluv = cieluv.CIELUV(whitepoint=whitepoint)
        return

    def from_xyz100(self, xyz):
        L, u, v = self.cieluv.from_xyz100(xyz)
        C = numpy.sqrt(u**2 + v**2)
        h = numpy.mod(numpy.arctan2(v, u), 2*numpy.pi) / numpy.pi * 180
        return numpy.array([L, C, h])

    def to_xyz100(self, lch):
        L, C, h = lch
        h_ = h * numpy.pi / 180
        luv = numpy.array([L, C * numpy.cos(h_), C * numpy.sin(h_)])
        return self.cieluv.to_xyz100(luv)
