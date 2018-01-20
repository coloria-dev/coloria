# -*- coding: utf-8 -*-
#
from __future__ import division

from .srgb_linear import SrgbLinear
from .srgb1 import SRGB1


def srgb_gamut(filename='srgb-xyz.vtu', n=50):
    import meshio
    import meshzoo
    points, cells = meshzoo.cube(nx=n, ny=n, nz=n)
    pts = SrgbLinear().to_xyz(points.T).T
    meshio.write(
        filename,
        pts, {'tetra': cells},
        point_data={'srgb1': SRGB1().from_srgb_linear(points)}
        )
    return
