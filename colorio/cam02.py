# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .illuminants import white_point, d65
from .ciecam02 import CIECAM02
from .srgb_linear import SrgbLinear
from .srgb1 import SRGB1


class CAM02_UCS(object):
    def __init__(self, c, Y_b, L_A, whitepoint=white_point(d65())):
        self.ciecam02 = CIECAM02(c, Y_b, L_A, whitepoint)
        self.K_L = 1.0
        self.c_1 = 0.007
        self.c_2 = 0.0228
        return

    def from_xyz(self, xyz):
        J, _, _, h, M, _, _ = self.ciecam02.from_xyz(xyz)
        J_ = (1+100*self.c_1)*J / (1 + self.c_1*J)
        M_ = 1/self.c_2 * numpy.log(1 + self.c_2*M)
        return numpy.array([J_, M_*numpy.cos(h), M_*numpy.sin(h)])

    def to_xyz(self, jab):
        J_, a, b = jab

        J = J_ / (1 - (J_-100)*self.c_1)

        h = numpy.arctan2(b, a)
        M_ = numpy.sqrt(a**2 + b**2)
        M = (numpy.exp(M_ * self.c_2) - 1) / self.c_2

        return self.ciecam02.to_xyz(numpy.array([J, M, h]), 'JMh')

    def srgb_gamut(self, filename='srgb-cam02ucs.vtu', n=50):
        import meshio
        import meshzoo
        points, cells = meshzoo.cube(nx=n, ny=n, nz=n)
        pts = self.from_xyz(SrgbLinear().to_xyz(points.T)).T
        rgb = SRGB1().from_srgb_linear(points)
        meshio.write(
            filename,
            pts, {'tetra': cells},
            point_data={'srgb': rgb}
            )
        return
