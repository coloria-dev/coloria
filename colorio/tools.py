# -*- coding: utf-8 -*-
#
from __future__ import division

import os

import matplotlib
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy
from scipy.optimize import leastsq
from scipy.spatial import ConvexHull
import yaml

from .illuminants import (
    spectrum_to_xyz100, planckian_radiator, whitepoints_cie1931
    )
from . import observers
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


def _plot_monochromatic(observer, xy_to_2d):
    # draw outline of monochromatic spectra
    lmbda = 1.0e-9 * numpy.arange(380, 701)
    values = []
    # TODO vectorize (see <https://github.com/numpy/numpy/issues/10439>)
    xyy = XYY()
    for k, _ in enumerate(lmbda):
        data = numpy.zeros(len(lmbda))
        data[k] = 1.0
        values.append(xy_to_2d(
            xyy.from_xyz100(spectrum_to_xyz100((lmbda, data), observer))[:2]
            ))
    values = numpy.array(values)
    # fill horseshoe area
    plt.fill(values[:, 0], values[:, 1], color=[0.8, 0.8, 0.8], zorder=0)
    # plot horseshoe outline
    plt.plot(values[:, 0], values[:, 1], '-k', label='monochromatic light')
    return


def _plot_rgb_triangle(xy_to_2d, bright=True):
    # plot sRGB triangle
    # discretization points
    n = 50

    # Get all RGB values that sum up to 1.
    rgb_linear = numpy.array(partition(3, n)).T / n
    if bright:
        # For the x-y-diagram, it doesn't matter if the values are scaled in
        # any way. After all, the tranlation to XYZ is linear, and then to xyY
        # it's (X/(X+Y+Z), Y/(X+Y+Z), Y), so the factor will only be present in
        # the last component which is discarded. To make the plot a bit
        # brighter, scale the colors up as much as possible.
        rgb_linear /= numpy.max(rgb_linear, axis=0)

    xyy = XYY()

    srgb_linear = SrgbLinear()
    xyz = srgb_linear.to_xyz100(rgb_linear)
    xyy_vals = xy_to_2d(xyy.from_xyz100(xyz)[:2])

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


def _plot_planckian_locus(observer, xy_to_2d):
    # plot planckian locus
    xyy = XYY()
    values = []
    for temp in numpy.arange(1000, 20001, 100):
        xyy_vals = xy_to_2d(xyy.from_xyz100(
            spectrum_to_xyz100(planckian_radiator(temp), observer)
            ))
        values.append(xyy_vals)
    values = numpy.array(values)
    plt.plot(values[:, 0], values[:, 1], ':k', label='Planckian locus')
    return


def plot_flat_gamut(
        xy_to_2d=lambda xy: xy,
        axes_labels=['x', 'y'],
        plot_rgb_triangle=True, plot_planckian_locus=True,
        ):
    '''Show the (u', v') gamut from CIELUV. There exists a chroma gamut for
    this color model as lines in XYZ are transformed to lines in CIELUV, hence
    CIELUV has a natural decomposition into lightness and chroma components.
    Also, the flat gamut is the same for every lightness value; this is _not_
    the case for most other color spaces (e.g., CIELAB).
    '''
    observer = observers.cie_1931_2()
    # observer = observers.cie_1964_10()

    _plot_monochromatic(observer, xy_to_2d)
    # plt.grid()

    if plot_rgb_triangle:
        _plot_rgb_triangle(xy_to_2d)
    if plot_planckian_locus:
        _plot_planckian_locus(observer, xy_to_2d)

    plt.gca().set_aspect('equal')
    # plt.legend()
    plt.xlabel('u\'')
    plt.ylabel('v\'')
    return


def show_flat_gamut(*args, **kwargs):
    plot_flat_gamut(*args, **kwargs)
    plt.show()
    return


def show_ebner_fairchild(colorspace):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, 'data/ebner_fairchild.yaml')) as f:
        data = yaml.safe_load(f)

    wp = colorspace.from_xyz100(numpy.array(data['white point']))[1:]

    d = [
        numpy.column_stack([dat['reference xyz'], numpy.array(dat['same']).T])
        for dat in data['data']
        ]

    _show_color_constancy_data(d, wp, colorspace)
    return


def show_hung_berns(colorspace):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, 'data/hung-berns/table3.yaml')) as f:
        data = yaml.safe_load(f)

    wp = colorspace.from_xyz100(numpy.array(whitepoints_cie1931['C']))[1:]

    d = [numpy.array(list(color.values())).T for color in data.values()]

    _show_color_constancy_data(d, wp, colorspace)
    return


def show_xiao(colorspace):
    dir_path = os.path.dirname(os.path.realpath(__file__))

    filenames = [
        'unique_blue.yaml', 'unique_green.yaml', 'unique_red.yaml',
        'unique_yellow.yaml'
        ]

    data = []
    for filename in filenames:
        with open(os.path.join(dir_path, 'data', 'xiao', filename)) as f:
            dat = numpy.array(yaml.safe_load(f))
        # average over all observers and sessions
        data.append(
            numpy.sum(dat, axis=(0, 1)) / numpy.prod(dat.shape[:2])
            )

    data = numpy.array(data)

    # Use Xiao's 'neutral gray' as white point.
    with open(os.path.join(dir_path, 'data/xiao/neutral_gray.yaml')) as f:
        ng_data = numpy.array(yaml.safe_load(f))

    ng = numpy.sum(ng_data, axis=0) / numpy.prod(ng_data.shape[:1])
    ng_cs = colorspace.from_xyz100(ng)[1:]

    data = numpy.moveaxis(data, 1, 2)

    _show_color_constancy_data(data, ng_cs, colorspace)
    return


def _show_color_constancy_data(data, wp, colorspace):
    srgb = SrgbLinear()
    for xyz in data:
        d = colorspace.from_xyz100(xyz)[1:]

        # Find best fit line through all points
        def f(theta, D=d):
            return (
                + numpy.sin(theta) * (D[0] - wp[0])
                + numpy.cos(theta) * (D[1] - wp[1])
                )

        def jac(theta, D=d):
            return (
                + numpy.cos(theta) * (D[0] - wp[0])
                - numpy.sin(theta) * (D[1] - wp[1])
                )

        # out = least_squares(f, 0.0)
        # out = leastsq(f, 0.0, full_output=True)

        out, _ = leastsq(f, 0.0, Dfun=jac)
        # We have to take the first element here, see
        # <https://github.com/scipy/scipy/issues/8532>.
        theta = out[0]

        # Plot it from wp to the outmost point
        length = numpy.sqrt(numpy.max(
            numpy.einsum('ij,ij->i', (d.T-wp), (d.T-wp))
            ))
        # The solution theta can be rotated by pi and still give the same
        # result. Find out on which side all the points are sitting and plot
        # the line accordingly.
        ex = length * numpy.array([numpy.cos(theta), -numpy.sin(theta)])

        end_point = wp + ex
        ep_d = numpy.linalg.norm(end_point - d[:, -1])
        ep_wp = numpy.linalg.norm(end_point - wp)
        if ep_d > ep_wp:
            end_point = wp - ex
        plt.plot(
            [wp[0], end_point[0]], [wp[1], end_point[1]], '-', color='0.5'
            )

        # Deliberatly only handle the two last components, e.g., a* b* from
        # L*a*b*. They typically indicate the chroma.
        for dd, rgb in zip(d.T, srgb.from_xyz100(xyz).T):
            is_legal_srgb = numpy.all(rgb >= 0) and numpy.all(rgb <= 1)
            col = srgb.to_srgb1(rgb) if is_legal_srgb else 'white'
            ecol = srgb.to_srgb1(rgb) if is_legal_srgb else 'black'
            plt.plot(dd[0], dd[1], 'o', color=col, markeredgecolor=ecol)

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
                 plot_standard_deviations=False,
                 xy_to_2d=lambda xy: xy
                 ):
    '''See <https://en.wikipedia.org/wiki/MacAdam_ellipse>,
    <https://doi.org/10.1364%2FJOSA.32.000247>.
    '''
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, 'data/macadam1942/table3.yaml')) as f:
        data = yaml.safe_load(f)

    plot_flat_gamut(plot_planckian_locus=False, xy_to_2d=xy_to_2d)
    # plt.grid(zorder=0)
    ax = plt.gca()

    # if plot_filter_positions:
    #     with open(os.path.join(dir_path, 'data/macadam1942/table1.yaml')) as f:
    #         filters_xyz = yaml.safe_load(f)
    #     filters_xyz = {
    #         key: 100 * numpy.array(value) for key, value in filters_xyz.items()
    #         }
    #     for key, xyz in filters_xyz.items():
    #         x, y = xyz100_to_2d(xyz)
    #         plt.plot(x, y, 'xk')
    #         ax.annotate(key, (x, y))

    for datak in data:
        # collect ellipse points
        _, _, _, _, delta_y_delta_x, delta_s = numpy.array(datak['data']).T

        center = xy_to_2d([datak['x'], datak['y']])

        X = (xy_to_2d(
            (center +
            (numpy.array([numpy.ones(delta_y_delta_x.shape[0]), delta_y_delta_x])
            / numpy.sqrt(1 + delta_y_delta_x**2) * delta_s).T).T
            ).T - xy_to_2d(center)).T

        if X.shape[1] < 2:
            continue

        # # Curve fit an ellipse with the data
        # (a, b, theta), _ = curve_fit(
        #     # dont' divide by a**2, b**2 here to avoid numerical difficulties
        #     # when optimizing
        #     lambda X, a, b, theta: (
        #         + a**2 * (X[0] * numpy.cos(theta) + X[1] * numpy.sin(theta))**2
        #         + b**2 * (X[0] * numpy.sin(theta) - X[1] * numpy.cos(theta))**2
        #         ),
        #     X, numpy.ones(X.shape[1])
        #     )

        def f_ellipse(a_b_theta, x=X):
            a, b, theta = a_b_theta
            return (
                + a**2 * (x[0] * numpy.cos(theta) + x[1] * numpy.sin(theta))**2
                + b**2 * (x[0] * numpy.sin(theta) - x[1] * numpy.cos(theta))**2
                - 1.0
                )

        def jac(a_b_theta, x=X):
            a, b, theta = a_b_theta
            return numpy.array([
                + 2*a * (x[0] * numpy.cos(theta) + x[1] * numpy.sin(theta))**2,
                + 2*b * (x[0] * numpy.sin(theta) - x[1] * numpy.cos(theta))**2,
                + a**2 * 2*(x[0] * numpy.cos(theta) + x[1] * numpy.sin(theta))
                * (-x[0] * numpy.sin(theta) + x[1] * numpy.cos(theta))
                + b**2 * 2*(x[0] * numpy.sin(theta) - x[1] * numpy.cos(theta))
                * (x[0] * numpy.cos(theta) + x[1] * numpy.sin(theta)),
                ]).T

        # out = leastsq(f_ellipse, [1.0, 1.0, 0.0], Dfun=jac, full_output=True)
        # print(out)
        (a, b, theta), _ = leastsq(f_ellipse, [1.0, 1.0, 0.0], Dfun=jac)
        # (a, b, theta), _, infodict, msg, ierr = \
        #     leastsq(f, [1.0, 1.0, 0.0], full_output=True, Dfun=jac)
        # print(infodict['nfev'])

        if plot_standard_deviations:
            # Original unscaled ellipses with data points
            plt.plot(*(X.T + center).T, '.', color='k')
            plt.plot(*center, 'x', color='k')
            e = Ellipse(
                xy=center,
                width=2/a, height=2/b,
                angle=theta / numpy.pi * 180
                )
            ax.add_artist(e)
            e.set_alpha(0.5)
            e.set_facecolor('k')

        # plot the scaled ellipse
        e = Ellipse(
            xy=center,
            width=scaling*2/a, height=scaling*2/b,
            angle=theta / numpy.pi * 180
            )
        ax.add_artist(e)
        e.set_alpha(0.5)
        e.set_facecolor('k')

    plt.show()
    return


def show_luo_rigg(scaling=1):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, 'data/luo-rigg/luo-rigg.yaml')) as f:
        data = yaml.safe_load(f)

    plot_xy_gamut(plot_planckian_locus=False, plot_rgb_triangle=False)
    ax = plt.gca()

    for _, data_set in data.items():
        # The set factor is the mean of the R values
        # set_factor = (
        #     numpy.sum(numpy.array(list(data_set.values()))[:, -1])
        #     / len(data_set)
        #     )
        for _, dat in data_set.items():
            x, y, Y, a, ab, theta, _ = dat
            a /= 1.0e4
            a *= (Y/30)**0.2

            b = a / ab

            a *= scaling  # * R / set_factor
            b *= scaling  # * R / set_factor

            # plot the ellipse
            e = Ellipse(xy=[x, y], width=2*a, height=2*b, angle=theta)
            ax.add_artist(e)
            e.set_alpha(0.5)
            e.set_facecolor('black')

    plt.show()
    return


def show_straights(cs):
    from mpl_toolkits.mplot3d import Axes3D
    # Some straight lines in XYZ
    t = numpy.linspace(0.0, 1.0, 101)
    n = 10

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # ax.set_aspect('equal')

    for k in range(n):
        s1 = numpy.random.rand(3)
        s1 /= numpy.linalg.norm(s1)
        s1 *= 100
        line = numpy.outer(s1, t)
        cs_line = cs.from_xyz100(line)
        # ax.plot(
        #     [cs_line[0][0], cs_line[0][-1]],
        #     [cs_line[1][0], cs_line[1][-1]],
        #     [cs_line[2][0], cs_line[2][-1]],
        #     color='0.5'
        #     )
        ax.plot(*cs_line)

    ax.set_xlabel(cs.labels[0])
    ax.set_ylabel(cs.labels[1])
    ax.set_zlabel(cs.labels[2])

    plt.show()
    return
