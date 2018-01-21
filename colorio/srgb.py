# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy


class SrgbLinear(object):
    def __init__(self):
        self.M = numpy.array([
            [+3.2406255, -1.537208, -0.4986286],
            [-0.9689307, +1.8757561, +0.0415175],
            [+0.0557101, -0.2040211, +1.0569959],
            ])
        return

    def from_xyz(self, xyz):
        # https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
        # http://www.color.org/srgb.pdf
        return numpy.dot(self.M, xyz)

    def to_xyz(self, srgb1_linear):
        return numpy.linalg.solve(self.M, srgb1_linear)

    def srgb_gamut(self, filename='srgb.vtu', n=50):
        show_gamut(filename, self.from_xyz, n=n)
        return


class SRGB1(object):
    def __init__(self):
        self.a = 0.055
        return

    def from_srgb_linear(self, srgb_linear):
        is_smaller = srgb_linear <= 0.0031308

        srgb = numpy.array(srgb_linear, dtype=float)
        srgb[is_smaller] *= 12.92
        srgb[~is_smaller] = (1+self.a) * srgb[~is_smaller]**(1/2.4) - self.a
        return srgb

    def to_srgb_linear(self, srgb1):
        srgb_linear = numpy.array(srgb1, dtype=float)

        # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
        is_smaller = srgb_linear <= 12.92 * 0.0031308  # 0.040449936

        srgb_linear[is_smaller] /= 12.92
        srgb_linear[~is_smaller] = (
            (srgb_linear[~is_smaller] + self.a) / (1+self.a)
            )**2.4
        return srgb_linear


def show_gamut(filename, from_xyz, n=50, cut_000=False):
    import meshio
    import meshzoo
    points, cells = meshzoo.cube(nx=n, ny=n, nz=n)

    if cut_000:
        # cut off [0, 0, 0] to avoid division by 0 in the xyz conversion
        points = points[1:]
        cells = cells[~numpy.any(cells == 0, axis=1)]
        cells -= 1

    pts = from_xyz(SrgbLinear().to_xyz(points.T)).T
    rgb = SRGB1().from_srgb_linear(points)
    meshio.write(
        filename,
        pts, {'tetra': cells},
        point_data={'srgb': rgb}
        )
    return
