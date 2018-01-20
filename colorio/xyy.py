# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from . import srgb_linear
from . import srgb1


def from_xyz(xyz):
    sum_xyz = numpy.sum(xyz, axis=0)
    x, y, _ = xyz
    return numpy.array([x/sum_xyz, y/sum_xyz, y])


def to_xyz(xyy):
    x, y, Y = xyy
    return numpy.array([Y/y*x, Y, Y/y * (1-x-y)])


def srgb_gamut(filename='srgb-xyy.vtu', n=50):
    import meshio
    import meshzoo
    points, cells = meshzoo.cube(nx=n, ny=n, nz=n)
    # cut off [0, 0, 0] to avoid division by 0 in the xyz conversion
    points = points[1:]
    cells = cells[~numpy.any(cells == 0, axis=1)]
    cells -= 1

    xyz = srgb_linear.to_xyz(points.T)
    pts = from_xyz(xyz).T
    rgb = srgb1.from_srgb_linear(points)
    meshio.write(
        filename,
        pts, {'tetra': cells},
        point_data={'srgb1': rgb}
        )
    return
