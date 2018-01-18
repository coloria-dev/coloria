# -*- coding: utf-8 -*-
#
from __future__ import division

import matplotlib
import matplotlib.pyplot as plt
import numpy

from .conversions import (
    spectrum_to_xyz, xyz_to_xyy, srgb1_to_xyz, xyz_to_srgb1, xyy_to_xyz
    )
from .illuminants import planckian_radiator


def show_gamut_diagram(*args, **kwargs):
    plot_gamut_diagram(*args, **kwargs)
    plt.show()
    return


def _partition(boxes, balls):
    # <https://stackoverflow.com/a/36748940/353337>
    def rec(boxes, balls, parent=tuple()):
        if boxes > 1:
            for i in range(balls + 1):
                for x in rec(boxes - 1, i, parent + (balls - i,)):
                    yield x
        else:
            yield parent + (balls,)

    return list(rec(boxes, balls))


def _plot_rgb_triangle():
    # plot sRGB triangle
    # discretization points
    n = 50
    corners = xyz_to_xyy(srgb1_to_xyz([[1, 0, 0], [0, 1, 0], [0, 0, 1]])).T

    bary = numpy.array(_partition(3, n)).T / n
    xyy = numpy.sum([
        numpy.outer(bary[k], corners[k]) for k in range(3)
        ], axis=0).T
    rgb = xyz_to_srgb1(xyy_to_xyz(xyy))
    # Some values can be slightly off (in the range of 1.0e-15)
    print(rgb-1.0)
    print(max(rgb.flatten()))
    assert numpy.all(rgb > -1.0e-14)
    # assert numpy.all(rgb-1.0 < 1.0e-2)
    rgb[rgb < 0] = 0.0
    rgb[rgb > 1] = 1.0

    # plt.plot(X[0], X[1], 'xk')

    # Unfortunately, one cannot yet use tripcolors with explicit RGB
    # specification (see
    # <https://github.com/matplotlib/matplotlib/issues/10265>). As a
    # workaround, associate range(n) data with the points and create a colormap
    # that associates the integer values with the respective RGBs.
    z = numpy.arange(xyy.shape[1])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            'gamut', rgb.T, N=rgb.shape[1]
            )

    triang = matplotlib.tri.Triangulation(xyy[0], xyy[1])
    plt.tripcolor(triang, z, shading='gouraud', cmap=cmap)
    return


def plot_gamut_diagram():
    # draw outline of monochromatic spectra
    lmbda = 1.0e-9 * numpy.arange(380, 701)
    values = []
    for k, wave_length in enumerate(lmbda):
        data = numpy.zeros(len(lmbda))
        data[k] = 1.0
        xyz = spectrum_to_xyz((lmbda, data))
        xyy = xyz_to_xyy(xyz)
        values.append(xyy[:2])
    values = numpy.array(values)
    # fill horseshoe area
    plt.fill(values[:, 0], values[:, 1], color=[0.8, 0.8, 0.8], zorder=0)
    # plot horseshoe outline
    plt.plot(values[:, 0], values[:, 1], '-k', label='monochromatic light')

    _plot_rgb_triangle()

    # plot planckian locus
    values = []
    for temp in [k*1000 for k in range(1, 11)]:
        lmbda, data = planckian_radiator(temp)
        xyz = xyz_to_xyy(spectrum_to_xyz((lmbda, data)))
        xyy = xyz_to_xyy(xyz)
        values.append(xyy[:2])
    values = numpy.array(values)
    plt.plot(values[:, 0], values[:, 1], ':k', label='Planckian locus')

    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.gca().set_aspect('equal')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    return
