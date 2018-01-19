# -*- coding: utf-8 -*-
#
from __future__ import division

import matplotlib
import matplotlib.pyplot as plt
import numpy

from .illuminants import spectrum_to_xyz
from . import srgb_linear
from . import srgb1
from . import xyy
from .illuminants import planckian_radiator


def show_gamut_diagram(*args, **kwargs):
    plot_gamut_diagram(*args, **kwargs)
    plt.show()
    return


def partition(boxes, balls):
    # <https://stackoverflow.com/a/36748940/353337>
    def rec(boxes, balls, parent=tuple()):
        if boxes > 1:
            for i in range(balls + 1):
                for x in rec(boxes - 1, i, parent + (balls - i,)):
                    yield x
        else:
            yield parent + (balls,)

    return list(rec(boxes, balls))


def _plot_monochromatic():
    # draw outline of monochromatic spectra
    lmbda = 1.0e-9 * numpy.arange(380, 701)
    values = []
    # TODO vectorize
    for k, _ in enumerate(lmbda):
        data = numpy.zeros(len(lmbda))
        data[k] = 1.0
        values.append(xyy.from_xyz(spectrum_to_xyz((lmbda, data)))[:2])
    values = numpy.array(values)
    # fill horseshoe area
    plt.fill(values[:, 0], values[:, 1], color=[0.8, 0.8, 0.8], zorder=0)
    # plot horseshoe outline
    plt.plot(values[:, 0], values[:, 1], '-k', label='monochromatic light')
    return


def _plot_rgb_triangle():
    # plot sRGB triangle
    # discretization points
    n = 50

    # Get all RGB values that sum up to 1.
    rgb_linear = numpy.array(partition(3, n)).T / n
    xyz = srgb_linear.to_xyz(rgb_linear)
    xyy_vals = xyy.from_xyz(xyz)

    # For each point in the x-y-diagram that is also an SRGB value, there are
    # many different SRGB values. Essentially, from [x, y, Y], the last
    # components can be chosen at will, and many choices will lead to a valid
    # SRGB. One choice would be to take the Y that results in the brightest
    # SRGB. For simplicity, we rather just take in the input SRGB here.
    rgb = srgb1.from_srgb_linear(rgb_linear)

    # Unfortunately, one cannot yet use tripcolors with explicit RGB
    # specification (see
    # <https://github.com/matplotlib/matplotlib/issues/10265>). As a
    # workaround, associate range(n) data with the points and create a colormap
    # that associates the integer values with the respective RGBs.
    z = numpy.arange(xyy_vals.shape[1])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'gamut', rgb.T, N=len(rgb.T)
        )

    triang = matplotlib.tri.Triangulation(xyy_vals[0], xyy_vals[1])
    plt.tripcolor(triang, z, shading='gouraud', cmap=cmap)
    return


def _plot_planckian_locus():
    # plot planckian locus
    values = []
    for temp in numpy.arange(1000, 20001, 100):
        xyy_vals = xyy.from_xyz(spectrum_to_xyz(planckian_radiator(temp)))
        values.append(xyy_vals[:2])
    values = numpy.array(values)
    plt.plot(values[:, 0], values[:, 1], ':k', label='Planckian locus')
    return


def plot_gamut_diagram():
    _plot_monochromatic()
    _plot_rgb_triangle()
    _plot_planckian_locus()

    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.gca().set_aspect('equal')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    return
