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


# https://stackoverflow.com/a/41856550/353337
def colors_to_cmap(colors):
    colors = numpy.asarray(colors)
    if colors.shape[1] == 3:
        colors = numpy.hstack((colors, numpy.ones((len(colors), 1))))
    steps = (0.5 + numpy.asarray(range(len(colors)-1), dtype=numpy.float))/(len(colors) - 1)
    return matplotlib.colors.LinearSegmentedColormap(
        'auto_cmap',
        {clrname: (
           [(0, col[0], col[0])] +
           [(step, c0, c1) for (step, c0, c1) in zip(steps, col[:-1], col[1:])] +
           [(1, col[-1], col[-1])]
           )
         for (clridx, clrname) in enumerate(['red', 'green', 'blue', 'alpha'])
         for col in [colors[:, clridx]]
         },
        N=len(colors)
        )


def _plot_rgb_triangle():
    # plot sRGB triangle
    # discretization points
    n = 50
    corners = numpy.array([
        xyz_to_xyy(srgb1_to_xyz(rgb))
        for rgb in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        ])
    bary = numpy.array(_partition(3, n)).T / n
    X = numpy.sum([
        numpy.outer(bary[k], corners[k]) for k in range(3)
        ], axis=0).T
    x = numpy.array(X[0])
    y = numpy.array(X[1])
    xyy = numpy.array([x, y, 1-x-y])
    rgb = xyz_to_srgb1(xyy_to_xyz(xyy))
    # plot the points
    # plt.plot(X[0], X[1], 'xk')
    z = numpy.arange(len(x))
    cmap = colors_to_cmap(rgb.T)

    triang = matplotlib.tri.Triangulation(x, y)
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
    # fill area
    plt.fill(values[:, 0], values[:, 1], color=[0.8, 0.8, 0.8], zorder=0)
    # plot outline
    plt.plot(values[:, 0], values[:, 1], '-k', label='monochromatic light')

    _plot_rgb_triangle()

    # # plot planckian locus
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
