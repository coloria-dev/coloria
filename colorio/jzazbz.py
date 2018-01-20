# -*- coding: utf-8 -*-
#
'''
Muhammad Safdar, Guihua Cui, Youn Jin Kim, and Ming Ronnier Luo,
Perceptually uniform color space for image signals including high dynamic range
and wide gamut,
Optics Express Vol. 25, Issue 13, pp. 15131-15151 (2017),
<https://doi.org/10.1364/OE.25.015131>.
'''
from __future__ import division

import numpy

from .illuminants import white_point, d65
from . import srgb_linear
from . import srgb1

b = 1.15
g = 0.66
c1 = 3424 / 2**12
c2 = 2413 / 2**7
c3 = 2392 / 2**7
n = 2610 / 2**14
p = 1.7*2523 / 2**5
d = -0.56
d0 = 1.6295499532821566e-11

M1 = numpy.array([
    [0.41478972, 0.579999, 0.0146480],
    [-0.2015100, 1.120649, 0.0531008],
    [-0.0166008, 0.264800, 0.6684799],
    ])

M2 = numpy.array([
    [0.5, 0.5, 0],
    [3.524000, -4.066708, +0.542708],
    [0.199076, +1.096799, -1.295875],
    ])


def from_xyz(xyz, whitepoint=white_point(d65())):
    x, y, z = (xyz.T / whitepoint).T
    x_ = b*x - (b-1)*z
    y_ = g*y - (g-1)*x
    lms = numpy.dot(M1, numpy.array([x_, y_, z]))
    lms_ = ((c1 + c2*(lms/10000)**n) / (1 + c3*(lms/10000)**n))**p
    iz, az, bz = numpy.dot(M2, lms_)
    jz = (1+d) * iz / (1+d*iz) - d0
    return numpy.array([jz, az, bz])


def to_xyz(jzazbz, whitepoint=white_point(d65())):
    jz, az, bz = jzazbz
    iz = (jz + d0) / (1 + d - d*(jz+d0))
    lms_ = numpy.linalg.solve(M2, numpy.array([iz, az, bz]))
    lms = 10000 * ((c1 - lms_**(1/p)) / (c3*lms_**(1/p) - c2))**(1/n)
    x_, y_, z_ = numpy.linalg.solve(M1, lms)
    x = (x_ + (b-1)*z_) / b
    y = (y_ + (g-1)*x) / g
    return (numpy.array([x, y, z_]).T * whitepoint).T


def srgb_gamut(filename='srgb-jzazbz.vtu', n=50):
    import meshio
    import meshzoo
    points, cells = meshzoo.cube(nx=n, ny=n, nz=n)
    pts = from_xyz(srgb_linear.to_xyz(points.T)).T
    rgb = srgb1.from_srgb_linear(points)
    meshio.write(
        filename,
        pts, {'tetra': cells},
        point_data={'srgb': rgb}
        )
    return
