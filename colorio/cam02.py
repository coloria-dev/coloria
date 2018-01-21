# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .illuminants import white_point, d65
from .ciecam02 import CIECAM02
from . import srgb


def _JMh_to_Jab(J, M, h, c1, c2):
    # J, _, _, h, M, _, _ = self.ciecam02.from_xyz(xyz)
    J_ = (1+100*c1)*J / (1 + c1*J)
    M_ = 1/c2 * numpy.log(1 + c2*M)
    return numpy.array([J_, M_*numpy.cos(h), M_*numpy.sin(h)])


def _Jab_to_JMh(jab, c1, c2):
    J_, a, b = jab
    J = J_ / (1 - (J_-100)*c1)
    h = numpy.arctan2(b, a)
    M_ = numpy.sqrt(a**2 + b**2)
    M = (numpy.exp(M_ * c2) - 1) / c2
    return numpy.array([J, M, h])


class CAM02_LCD(object):
    def __init__(self, c, Y_b, L_A, whitepoint=white_point(d65())):
        self.ciecam02 = CIECAM02(c, Y_b, L_A, whitepoint)
        self.K_L = 0.77
        self.c1 = 0.007
        self.c2 = 0.0053
        return

    def from_xyz(self, xyz):
        J, _, _, h, M, _, _ = self.ciecam02.from_xyz(xyz)
        return _JMh_to_Jab(J, M, h, self.c1, self.c2)

    def to_xyz(self, jab):
        return self.ciecam02.to_xyz(_Jab_to_JMh(jab, self.c1, self.c2), 'JMh')

    def srgb_gamut(self, filename='srgb-cam02lcd.vtu', n=50):
        srgb.show_gamut(filename, self.from_xyz, n=n)
        return


class CAM02_SCD(object):
    def __init__(self, c, Y_b, L_A, whitepoint=white_point(d65())):
        self.ciecam02 = CIECAM02(c, Y_b, L_A, whitepoint)
        self.K_L = 1.24
        self.c1 = 0.007
        self.c2 = 0.0363
        return

    def from_xyz(self, xyz):
        J, _, _, h, M, _, _ = self.ciecam02.from_xyz(xyz)
        return _JMh_to_Jab(J, M, h, self.c1, self.c2)

    def to_xyz(self, jab):
        return self.ciecam02.to_xyz(_Jab_to_JMh(jab, self.c1, self.c2), 'JMh')

    def srgb_gamut(self, filename='srgb-cam02scd.vtu', n=50):
        srgb.show_gamut(filename, self.from_xyz, n=n)
        return


class CAM02_UCS(object):
    def __init__(self, c, Y_b, L_A, whitepoint=white_point(d65())):
        self.ciecam02 = CIECAM02(c, Y_b, L_A, whitepoint)
        self.K_L = 1.0
        self.c1 = 0.007
        self.c2 = 0.0228
        return

    def from_xyz(self, xyz):
        J, _, _, h, M, _, _ = self.ciecam02.from_xyz(xyz)
        return _JMh_to_Jab(J, M, h, self.c1, self.c2)

    def to_xyz(self, jab):
        return self.ciecam02.to_xyz(_Jab_to_JMh(jab, self.c1, self.c2), 'JMh')

    def srgb_gamut(self, filename='srgb-cam02ucs.vtu', n=50):
        srgb.show_gamut(filename, self.from_xyz, n=n)
        return
