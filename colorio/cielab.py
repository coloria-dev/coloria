# -*- coding: utf-8 -*-
#
from __future__ import division

import matplotlib.pyplot as plt
import numpy

from .illuminants import spectrum_to_xyz, white_point, d65


def from_xyz(xyz, whitepoint):
    def f(t):
        delta = 6.0/29.0
        out = numpy.array(t, dtype=float)
        is_greater = out > delta**3
        is_smaller = numpy.logical_not(is_greater)
        out[is_greater] = numpy.cbrt(out[is_greater])
        out[is_smaller] = out[is_smaller]/3/delta**2 + 4.0/29.0
        return out

    fxyz = f(xyz / whitepoint)
    return numpy.array([
        116 * fxyz[1] - 16,
        500 * (fxyz[0] - fxyz[1]),
        200 * (fxyz[1] - fxyz[2]),
        ])


def to_xyz(cielab, whitepoint):
    def f1(t):
        delta = 6.0/29.0
        out = numpy.array(t, dtype=float)
        is_greater = out > delta
        is_smaller = numpy.logical_not(is_greater)
        out[is_greater] = out[is_greater]**3
        out[is_smaller] = 3*delta**2 + (out[is_smaller] - 4.0/29.0)

    return whitepoint * f1(numpy.array([
        (cielab[0]+16)/116 + cielab[1]/500,
        (cielab[0]+16)/116,
        (cielab[0]+16)/116 - cielab[2]/500,
        ]))


def show_luminance_level(*args, **kwargs):
    plot_luminance_level(*args, **kwargs)
    plt.show()
    return


def plot_luminance_level(L):
    _plot_horseshoe(L)

    plt.title('L*={}'.format(L))
    plt.gca().set_aspect('equal')
    plt.legend()
    plt.xlabel('a*')
    plt.ylabel('b*')
    return


def _plot_rgb(L):
    return


def _plot_horseshoe(L):
    wp = 100 * white_point(d65())

    # The target Y is chosen such that the desired luminance is achieved.
    delta = 6.0 / 29.0
    y_target = (
        ((L+16)/116)**3 * wp[1] if L > 8 else
        3*L*delta**2*wp[1]
        )

    # draw outline of monochromatic spectra
    lmbda = 1.0e-9 * numpy.arange(380, 701)

    # TODO vectorize
    values = []
    for k, _ in enumerate(lmbda):
        data = numpy.zeros(len(lmbda))
        data[k] = 1.0
        xyz = spectrum_to_xyz((lmbda, data))
        # normalize to get the target luminance
        xyz *= y_target / xyz[1]
        vals = from_xyz(xyz, wp)
        assert abs(vals[0] - L) < 1.0e-13
        values.append(vals[1:])

    values = numpy.array(values)

    plt.fill(values[:, 0], values[:, 1], color=[0.8, 0.8, 0.8], zorder=0)
    plt.plot(values[:, 0], values[:, 1], '-k', label='monochromatic light')
    return
