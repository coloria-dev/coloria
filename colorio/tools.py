# -*- coding: utf-8 -*-
#
from __future__ import division

import os

import matplotlib
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy
from scipy.optimize import curve_fit
from scipy.spatial import ConvexHull
import yaml

from .illuminants import (
    spectrum_to_xyz100, planckian_radiator, whitepoints_cie1931
    )
from .rec2020 import Rec2020
from .srgb import SrgbLinear
from .xyy import XYY


def delta(a, b):
    '''Computes the distances between two colors or color sets. The shape of
    `a` and `b` must be equal.
    '''
    diff = a - b
    return numpy.einsum('i...,i...->...', diff, diff)


def show_visible_gamut(colorspace, observer, illuminant, filename,
                       cut_000=False):
    import meshio

    # The XYZ gamut is actually defined by an arbitrarily chosen maximum
    # intensity (here: 1). Then, all block spectra with this intensity are
    # mapped into XYZ space; they form the outer hull.
    lmbda, illu = illuminant
    values = []

    data = numpy.zeros(len(lmbda))
    values.append(
        spectrum_to_xyz100((lmbda, illu*data), observer=observer)
        )
    for width in range(1, len(lmbda)):
        data = numpy.zeros(len(lmbda))
        data[:width] = 1.0
        for _, _ in enumerate(lmbda):
            values.append(
                spectrum_to_xyz100((lmbda, illu*data), observer=observer)
                )
            data = numpy.roll(data, shift=1)
    data = numpy.ones(len(lmbda))
    values.append(
        spectrum_to_xyz100((lmbda, illu*data), observer=observer)
        )

    # scale the values such that the Y-coordinate of the white point has value
    # 100.
    values = numpy.array(values)
    values *= 100 / values[-1][1]

    cells = ConvexHull(values).simplices

    if cut_000:
        values = values[1:]
        cells = cells[~numpy.any(cells == 0, axis=1)]
        cells -= 1

    pts = colorspace.from_xyz100(values.T).T

    meshio.write(filename, pts, cells={'triangle': cells})
    return


def show_srgb_gamut(colorspace, filename, n=50, cut_000=False):
    import meshio
    import meshzoo
    points, cells = meshzoo.cube(nx=n, ny=n, nz=n)

    if cut_000:
        # cut off [0, 0, 0] to avoid division by 0 in the xyz conversion
        points = points[1:]
        cells = cells[~numpy.any(cells == 0, axis=1)]
        cells -= 1

    srgb_linear = SrgbLinear()
    pts = colorspace.from_xyz100(srgb_linear.to_xyz100(points.T)).T
    rgb = srgb_linear.to_srgb1(points)
    meshio.write(
        filename,
        pts, {'tetra': cells},
        point_data={'srgb': rgb}
        )
    return


def show_hdr_gamut(colorspace, filename, n=50, cut_000=False):
    import meshio
    import meshzoo
    points, cells = meshzoo.cube(nx=n, ny=n, nz=n)

    if cut_000:
        # cut off [0, 0, 0] to avoid division by 0 in the xyz conversion
        points = points[1:]
        cells = cells[~numpy.any(cells == 0, axis=1)]
        cells -= 1

    cs = Rec2020()
    pts = colorspace.from_xyz100(cs.to_xyz100(points.T)).T
    rgb = cs.to_gamma(points)
    meshio.write(
        filename,
        pts, {'tetra': cells},
        point_data={'rec2020-rgb': rgb}
        )
    return


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
    # TODO vectorize (see <https://github.com/numpy/numpy/issues/10439>)
    for k, _ in enumerate(lmbda):
        data = numpy.zeros(len(lmbda))
        data[k] = 1.0
        values.append(XYY().from_xyz100(spectrum_to_xyz100((lmbda, data)))[:2])
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
    # For the x-y-diagram, it doesn't matter if the values are scaled in any
    # way. After all, the tranlation to XYZ is linear, and then to xyY it's
    # (X/(X+Y+Z), Y/(X+Y+Z), Y), so the factor will only be present in the last
    # component which is discarded. To make the plot a bit brighter, scale the
    # colors up as much as possible.
    rgb_linear /= numpy.max(rgb_linear, axis=0)

    srgb_linear = SrgbLinear()
    xyz = srgb_linear.to_xyz100(rgb_linear)
    xyy_vals = XYY().from_xyz100(xyz)

    # Unfortunately, one cannot use tripcolors with explicit RGB specification
    # (see <https://github.com/matplotlib/matplotlib/issues/10265>). As a
    # workaround, associate range(n) data with the points and create a colormap
    # that associates the integer values with the respective RGBs.
    z = numpy.arange(xyy_vals.shape[1])
    rgb = srgb_linear.to_srgb1(rgb_linear)
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
        xyy_vals = XYY().from_xyz100(
            spectrum_to_xyz100(planckian_radiator(temp))
            )
        values.append(xyy_vals[:2])
    values = numpy.array(values)
    plt.plot(values[:, 0], values[:, 1], ':k', label='Planckian locus')
    return


def plot_gamut_diagram(plot_rgb_triangle=True, plot_planckian_locus=True):
    plt.plot([0.0, 1.0], [1.0, 0.0], color='0.8')

    _plot_monochromatic()
    if plot_rgb_triangle:
        _plot_rgb_triangle()
    if plot_planckian_locus:
        _plot_planckian_locus()

    plt.xlim(xmin=0, xmax=0.8)
    plt.ylim(ymin=0, ymax=0.9)

    plt.gca().set_aspect('equal')
    # plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    return


def show_ebner_fairchild(colorspace):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, 'data/ebner_fairchild.yaml')) as f:
        data = yaml.safe_load(f)

    # show white point
    wp = colorspace.from_xyz100(numpy.array(data['white point']))
    plt.plot(wp[1], wp[2], '.k')

    srgb = SrgbLinear()
    for item in data['data']:
        rgb = srgb.to_srgb1(srgb.from_xyz100(item['reference xyz']))
        xyz = numpy.array(item['same']).T
        # The points are sorted by the first components d[0] (typically
        # luminance).
        d = colorspace.from_xyz100(xyz)

        # Deliberatly only handle the two last components, e.g., a* b* from
        # L*a*b*. They typically indicate the chroma.
        plt.plot(d[1], d[2], '-', color='0.5')
        for dd, xyz_ in zip(d.T, xyz.T):
            rgb = srgb.from_xyz100(xyz_)
            is_legal_srgb = numpy.all(rgb >= 0) and numpy.all(rgb <= 1)
            col = srgb.to_srgb1(rgb) if is_legal_srgb else 'white'
            ecol = srgb.to_srgb1(rgb) if is_legal_srgb else 'k'
            plt.plot(dd[1], dd[2], 'o', color=col, markeredgecolor=ecol)

    plt.axis('equal')
    plt.show()
    return


def show_hung_berns(colorspace):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, 'data/hung-berns/table3.yaml')) as f:
        data = yaml.safe_load(f)

    # show white point
    d = colorspace.from_xyz100(numpy.array(whitepoints_cie1931['C']))
    plt.plot(d[1], d[2], '.k')

    srgb = SrgbLinear()
    for color_name in data.keys():
        dat = data[color_name]
        xyz = numpy.array(list(dat.values())).T
        d = colorspace.from_xyz100(xyz)

        # Deliberatly only handle the two last components, e.g., a* b* from
        # L*a*b*. They typically indicate the chroma.
        # Plot the lines in black first, then the individual points.
        plt.plot(d[1], d[2], '-', color='k')
        for dd, rgb in zip(d.T, srgb.from_xyz100(xyz).T):
            is_legal_srgb = numpy.all(rgb >= 0) and numpy.all(rgb <= 1)
            col = srgb.to_srgb1(rgb) if is_legal_srgb else 'white'
            ecol = srgb.to_srgb1(rgb) if is_legal_srgb else 'black'
            plt.plot(dd[1], dd[2], 'o', color=col, markeredgecolor=ecol)

    plt.axis('equal')
    plt.show()
    return


def show_munsell(colorspace, V):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, 'data/munsell/real.yaml')) as f:
        data = yaml.safe_load(f)

    # https://stackoverflow.com/a/6473724/353337
    v = data['V']
    x = data['x']
    y = data['y']
    Y = data['Y']

    idx = numpy.array(v) == V
    x = numpy.array(x)
    y = numpy.array(y)
    Y = numpy.array(Y)
    xyy = numpy.array([x[idx], y[idx], Y[idx]])

    xyz = XYY().to_xyz100(xyy)
    vals = colorspace.from_xyz100(xyz)

    srgb = SrgbLinear()
    rgb = srgb.from_xyz100(xyz)
    is_legal_srgb = numpy.logical_and(
        numpy.all(rgb >= 0, axis=0), numpy.all(rgb <= 1, axis=0)
        )

    # plot the ones that cannot be represented in SRGB
    plt.plot(
        vals[1, ~is_legal_srgb], vals[2, ~is_legal_srgb],
        'o', color='white', markeredgecolor='black'
        )
    # plot the srgb dots
    for val, rgb_ in zip(vals[:, is_legal_srgb].T, rgb[:, is_legal_srgb].T):
        plt.plot(val[1], val[2], 'o', color=srgb.to_srgb1(rgb_))

    plt.axis('equal')
    plt.show()
    return


def show_macadam(scaling=1,
                 plot_filter_positions=False,
                 plot_standard_deviations=True):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, 'data/macadam1942/table3.yaml')) as f:
        data = yaml.safe_load(f)

    plot_gamut_diagram(plot_planckian_locus=False)
    ax = plt.gca()

    xyy = XYY()

    if plot_filter_positions:
        with open(os.path.join(dir_path, 'data/macadam1942/table1.yaml')) as f:
            filters_xyz = yaml.safe_load(f)
        filters_xyz = {
            key: 100 * numpy.array(value) for key, value in filters_xyz.items()
            }
        for key, xyz in filters_xyz.items():
            x, y = xyy.from_xyz100(xyz)[:2]
            plt.plot(x, y, 'xk')
            ax.annotate(key, (x, y))

    for datak in data:
        # collect ellipse points
        dat = numpy.array(datak['data']).T
        X = (
            numpy.array([numpy.ones(dat[4].shape[0]), dat[4]])
            / numpy.sqrt(1 + dat[4]**2) * dat[5]
            )

        if X.shape[1] < 2:
            continue

        # Curve fit an ellipse with the data
        (a, b, theta), _ = curve_fit(
            # dont' divide by a**2, b**2 here to avoid numerical difficulties
            # when optimizing
            lambda X, a, b, theta: (
                + a**2 * (X[0] * numpy.cos(theta) + X[1] * numpy.sin(theta))**2
                + b**2 * (X[0] * numpy.sin(theta) - X[1] * numpy.cos(theta))**2
                ),
            X, numpy.ones(X.shape[1])
            )

        # plot the ellipse
        e = Ellipse(
            xy=[datak['x'], datak['y']],
            width=scaling * 2/a, height=scaling * 2/b,
            angle=theta / numpy.pi * 180
            )
        ax.add_artist(e)
        e.set_alpha(0.5)
        e.set_facecolor('k')

        if plot_standard_deviations:
            plt.plot(
                [datak['x'] - scaling*X[0], datak['x'] + scaling*X[0]],
                [datak['y'] - scaling*X[1], datak['y'] + scaling*X[1]],
                'kx'
                )

    plt.show()
    return
