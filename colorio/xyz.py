# -*- coding: utf-8 -*-
#
from __future__ import division

from . import srgb_linear
from . import srgb1


def srgb_gamut(filename='srgb-xyz.vtu', n=50):
    import meshio
    import meshzoo
    points, cells = meshzoo.cube(nx=n, ny=n, nz=n)
    pts = srgb_linear.to_xyz(points.T).T
    meshio.write(
        filename,
        pts, {'tetra': cells},
        point_data={'srgb1': srgb1.from_srgb_linear(points)}
        )
    return
