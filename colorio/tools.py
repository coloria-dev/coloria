# -*- coding: utf-8 -*-
#
from __future__ import division

import matplotlib.pyplot as plt
import numpy

from .conversions import spectrum_to_xyz, xyz_to_xyy, srgb1_to_xyz
from .illuminants import planckian_radiator


def show_gamut_diagram(*args, **kwargs):
    plot_gamut_diagram(*args, **kwargs)
    plt.show()
    return


def plot_gamut_diagram():
    # draw outline of monochromatic spectra
    lmbda = numpy.arange(400, 701, 1)
    values = []
    for k, wave_length in enumerate(lmbda):
        data = numpy.zeros(len(lmbda))
        data[k] = 1.0
        xyz = spectrum_to_xyz((lmbda, data))
        xyy = xyz_to_xyy(xyz)
        values.append(xyy[:2])
    values = numpy.array(values)
    plt.plot(values[:, 0], values[:, 1], '-k', label='monochromatic light')

    # plot sRGB triangle
    for rgb in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
        xyy = xyz_to_xyy(srgb1_to_xyz(rgb))
        plt.plot(xyy[0], xyy[1], 'x', color=rgb)

    # plot planckian locus
    values = []
    for temp in [k*1000 for k in range(1, 11)]:
        lmbda, data = planckian_radiator(temp)
        xyz = xyz_to_xyy(spectrum_to_xyz((lmbda, data)))
        xyy = xyz_to_xyy(xyz)
        print(temp, xyz, xyy)
        values.append(xyy[:2])
    values = numpy.array(values)
    plt.plot(values[:, 0], values[:, 1], '-x', label='Planckian locus')

    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.gca().set_aspect('equal')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    return
