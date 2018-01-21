# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .illuminants import white_point, d65
from .ciecam02 import CIECAM02


class CAM02(object):
    def __init__(self, variant, c, Y_b, L_A, whitepoint=white_point(d65())):
        params = {
            'LCD': (0.77, 0.007, 0.0053),
            'SCD': (1.24, 0.007, 0.0363),
            'UCS': (1.00, 0.007, 0.0228),
            }
        self.K_L, self.c1, self.c2 = params[variant]
        self.ciecam02 = CIECAM02(c, Y_b, L_A, whitepoint)
        return

    def from_xyz(self, xyz):
        J, _, _, h, M, _, _ = self.ciecam02.from_xyz(xyz)
        J_ = (1+100*self.c1)*J / (1 + self.c1*J)
        M_ = 1/self.c2 * numpy.log(1 + self.c2*M)
        return numpy.array([J_, M_*numpy.cos(h), M_*numpy.sin(h)])

    def to_xyz(self, jab):
        J_, a, b = jab
        J = J_ / (1 - (J_-100)*self.c1)
        h = numpy.arctan2(b, a)
        M_ = numpy.sqrt(a**2 + b**2)
        M = (numpy.exp(M_ * self.c2) - 1) / self.c2
        return self.ciecam02.to_xyz(numpy.array([J, M, h]), 'JMh')
