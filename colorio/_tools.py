import os

import matplotlib
import matplotlib.pyplot as plt
import numpy
import yaml
from matplotlib.patches import Ellipse

from . import observers
from ._srgb import SrgbLinear
from .illuminants import planckian_radiator, spectrum_to_xyz100


def delta(a, b):
    """Computes the distances between two colors or color sets. The shape of
    `a` and `b` must be equal.
    """
    diff = a - b
    return numpy.einsum("i...,i...->...", diff, diff)


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


def _xyy_from_xyz100(xyz):
    sum_xyz = numpy.sum(xyz, axis=0)
    x = xyz[0]
    y = xyz[1]
    return numpy.array([x / sum_xyz, y / sum_xyz, y / 100])


def _plot_monochromatic(observer, xy_to_2d, fill_horseshoe=True):
    # draw outline of monochromatic spectra
    lmbda = 1.0e-9 * numpy.arange(380, 701)
    values = []
    # TODO vectorize (see <https://github.com/numpy/numpy/issues/10439>)
    for k, _ in enumerate(lmbda):
        data = numpy.zeros(len(lmbda))
        data[k] = 1.0
        values.append(_xyy_from_xyz100(spectrum_to_xyz100((lmbda, data), observer))[:2])
    values = numpy.array(values)

    # Add the values between the first and the last point of the horseshoe
    t = numpy.linspace(0.0, 1.0, 101)
    connect = xy_to_2d(numpy.outer(values[0], t) + numpy.outer(values[-1], 1 - t))
    values = xy_to_2d(values.T).T
    full = numpy.concatenate([values, connect.T])

    # fill horseshoe area
    if fill_horseshoe:
        plt.fill(*full.T, color=[0.8, 0.8, 0.8], zorder=0)
    # plot horseshoe outline
    plt.plot(
        values[:, 0],
        values[:, 1],
        "-k",
        # label="monochromatic light"
    )
    # plot dotted connector
    plt.plot(connect[0], connect[1], ":k")
    return


def _plot_rgb_triangle(xy_to_2d, bright=True):
    # plot sRGB triangle
    # discretization points
    n = 50

    # Get all RGB values that sum up to 1.
    rgb_linear = numpy.array(partition(3, n)).T / n
    if bright:
        # For the x-y-diagram, it doesn't matter if the values are scaled in any way.
        # After all, the tranlation to XYZ is linear, and then to xyY it's (X/(X+Y+Z),
        # Y/(X+Y+Z), Y), so the factor will only be present in the last component which
        # is discarded. To make the plot a bit brighter, scale the colors up as much as
        # possible.
        rgb_linear /= numpy.max(rgb_linear, axis=0)

    srgb_linear = SrgbLinear()
    xyz = srgb_linear.to_xyz100(rgb_linear)
    xyy_vals = xy_to_2d(_xyy_from_xyz100(xyz)[:2])

    # Unfortunately, one cannot use tripcolors with explicit RGB specification
    # (see <https://github.com/matplotlib/matplotlib/issues/10265>). As a
    # workaround, associate range(n) data with the points and create a colormap
    # that associates the integer values with the respective RGBs.
    z = numpy.arange(xyy_vals.shape[1])
    rgb = srgb_linear.to_srgb1(rgb_linear)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "gamut", rgb.T, N=len(rgb.T)
    )

    triang = matplotlib.tri.Triangulation(xyy_vals[0], xyy_vals[1])
    plt.tripcolor(triang, z, shading="gouraud", cmap=cmap)
    return


def _plot_planckian_locus(observer, xy_to_2d):
    # plot planckian locus
    values = []
    for temp in numpy.arange(1000, 20001, 100):
        xyy_vals = xy_to_2d(
            _xyy_from_xyz100(spectrum_to_xyz100(planckian_radiator(temp), observer))
        )
        values.append(xyy_vals)
    values = numpy.array(values)
    plt.plot(values[:, 0], values[:, 1], ":k", label="Planckian locus")
    return


def plot_flat_gamut(
    xy_to_2d=lambda xy: xy,
    axes_labels=("x", "y"),
    plot_rgb_triangle=True,
    fill_horseshoe=True,
    plot_planckian_locus=True,
):
    """Show a flat color gamut, by default xy.  There exists a chroma gamut for
    all color models which transform lines in XYZ to lines, and hence have a
    natural decomposition into lightness and chroma components.  Also, the flat
    gamut is the same for every lightness value. Examples for color models with
    this property are CIELUV and IPT, examples for color models without are
    CIELAB and CIECAM02.
    """
    observer = observers.cie_1931_2()
    # observer = observers.cie_1964_10()

    _plot_monochromatic(observer, xy_to_2d, fill_horseshoe=fill_horseshoe)
    # plt.grid()

    if plot_rgb_triangle:
        _plot_rgb_triangle(xy_to_2d)
    if plot_planckian_locus:
        _plot_planckian_locus(observer, xy_to_2d)

    plt.gca().set_aspect("equal")
    # plt.legend()
    plt.xlabel(axes_labels[0])
    plt.ylabel(axes_labels[1])
    return


def show_flat_gamut(*args, **kwargs):
    plot_flat_gamut(*args, **kwargs)
    plt.show()
    return


def get_munsell_data():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, "data/munsell/real.yaml")) as f:
        data = yaml.safe_load(f)

    h = numpy.array(data["h"])
    V = numpy.array(data["V"])
    C = numpy.array(data["C"])
    xyy = numpy.array([data["x"], data["y"], data["Y"]])

    return h, V, C, xyy


def show_macadam(*args, **kwargs):
    """See <https://en.wikipedia.org/wiki/MacAdam_ellipse>,
    <https://doi.org/10.1364%2FJOSA.32.000247>.
    """
    plot_macadam(*args, **kwargs)
    plt.show()
    return


def save_macadam(filename, *args, **kwargs):
    plt.figure()
    plot_macadam(*args, **kwargs)
    plt.savefig(filename, bbox_inches="tight", transparent=True)
    return


def plot_macadam(
    ellipse_scaling=10,
    plot_filter_positions=False,
    plot_standard_deviations=False,
    plot_rgb_triangle=True,
    mesh_resolution=1,
    xy_to_2d=lambda xy: xy,
    axes_labels=("x", "y"),
):
    """See <https://en.wikipedia.org/wiki/MacAdam_ellipse>,
    <https://doi.org/10.1364%2FJOSA.32.000247>.
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, "data/macadam1942/table3.yaml")) as f:
        data = yaml.safe_load(f)

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

    # collect the ellipse centers and offsets
    centers = []
    offsets = []
    for datak in data:
        # collect ellipse points
        _, _, _, _, delta_y_delta_x, delta_s = numpy.array(datak["data"]).T

        offset = (
            numpy.array([numpy.ones(delta_y_delta_x.shape[0]), delta_y_delta_x])
            / numpy.sqrt(1 + delta_y_delta_x ** 2)
            * delta_s
        )

        if offset.shape[1] < 2:
            continue

        centers.append([datak["x"], datak["y"]])
        offsets.append(numpy.column_stack([+offset, -offset]))

    centers = numpy.array(centers)

    _plot_ellipse_data(
        centers,
        offsets,
        ellipse_scaling=ellipse_scaling,
        xy_to_2d=xy_to_2d,
        mesh_resolution=mesh_resolution,
        plot_rgb_triangle=plot_rgb_triangle,
    )
    return


def show_luo_rigg(*args, **kwargs):
    plot_luo_rigg(*args, **kwargs)
    plt.show()
    return


def save_luo_rigg(filename, *args, **kwargs):
    plt.figure()
    plot_luo_rigg(*args, **kwargs)
    plt.savefig(filename, bbox_inches="tight", transparent=True)
    return


def plot_luo_rigg(
    plot_rgb_triangle=True,
    ellipse_scaling=1,
    xy_to_2d=lambda xy: xy,
    mesh_resolution=1,
    ellipse_color="k",
):
    # M. R. Luo, B. Rigg,
    # Chromaticity Discrimination Ellipses for Surface Colours,
    # Color Research and Application, Volume 11, Issue 1, Spring 1986, Pages 25-42,
    # <https://doi.org/10.1002/col.5080110107>.

    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, "data/luo-rigg/luo-rigg.yaml")) as f:
        data = yaml.safe_load(f)

    centers = []
    offsets = []

    # collect the ellipse centers and offsets
    # Use four offset points of each ellipse, one could take more
    alpha = 2 * numpy.pi * numpy.linspace(0.0, 1.0, 16, endpoint=False)
    pts = numpy.array([numpy.cos(alpha), numpy.sin(alpha)])
    for data_set in data.values():
        # The set factor is the mean of the R values
        # set_factor = sum([dat[-1] for dat in data_set.values()]) / len(data_set)

        for dat in data_set.values():
            x, y, Y, a, a_div_b, theta_deg, _ = dat
            theta = theta_deg * 2 * numpy.pi / 360
            a /= 1.0e4
            a *= (Y / 30) ** 0.2
            b = a / a_div_b

            # a *= R / set_factor
            # b *= R / set_factor

            # plot the ellipse
            centers.append([x, y])
            J = numpy.array(
                [
                    [+a * numpy.cos(theta), -b * numpy.sin(theta)],
                    [+a * numpy.sin(theta), +b * numpy.cos(theta)],
                ]
            )
            offsets.append(numpy.dot(J, pts))

    centers = numpy.array(centers)
    _plot_ellipse_data(
        centers,
        offsets,
        ellipse_scaling=ellipse_scaling,
        mesh_resolution=mesh_resolution,
        xy_to_2d=xy_to_2d,
        plot_rgb_triangle=plot_rgb_triangle,
        ellipse_color=ellipse_color,
    )
    plt.xlim(0.0)
    plt.ylim(0.0)
    return


def _plot_ellipse_data(
    centers,
    offsets,
    xy_to_2d=lambda xy: xy,
    axes_labels=("x", "y"),
    plot_rgb_triangle=False,
    ellipse_scaling=10,
    ellipse_color="k",
    mesh_resolution=None,
):
    import meshzoo

    plot_flat_gamut(
        plot_planckian_locus=False,
        xy_to_2d=xy_to_2d,
        axes_labels=axes_labels,
        plot_rgb_triangle=plot_rgb_triangle,
        fill_horseshoe=mesh_resolution is None,
    )

    if mesh_resolution is not None:
        # dir_path = os.path.dirname(os.path.realpath(__file__))
        # with open(os.path.join(dir_path, 'data/gamut_triangulation.yaml')) as f:
        #     data = yaml.safe_load(f)
        # points = numpy.array(data['points'])
        # cells = numpy.array(data['cells'])

        points, cells = meshzoo.triangle(
            mesh_resolution,
            corners=numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        )
        points = points[:, :2]

        edges, _ = meshzoo.create_edges(cells)
        pts = xy_to_2d(points.T).T
        lines = pts[edges].T
        plt.plot(*lines, color="0.8", zorder=0)

    _plot_ellipses(centers, offsets, xy_to_2d, ellipse_scaling, facecolor=ellipse_color)
    return


def _plot_ellipses(
    centers, offsets, xy_to_2d, ellipse_scaling, alpha=0.5, facecolor="k", label=None
):
    from scipy.optimize import leastsq

    ax = plt.gca()

    for center, offset in zip(centers, offsets):
        # If xy_to_2d was linear, we would only need one of center+-offset
        tcenter = xy_to_2d(center)
        X = numpy.column_stack(
            [xy_to_2d((center + offset.T).T), xy_to_2d((center - offset.T).T)]
        )
        X = (X.T - tcenter).T

        def f_ellipse(a_b_theta, x=X):
            a, b, theta = a_b_theta
            sin_t = numpy.sin(theta)
            cos_t = numpy.cos(theta)
            return (
                +(a ** 2) * (x[0] * cos_t + x[1] * sin_t) ** 2
                + b ** 2 * (x[0] * sin_t - x[1] * cos_t) ** 2
                - 1.0
            )

        def jac(a_b_theta, x=X):
            a, b, theta = a_b_theta
            sin_t = numpy.sin(theta)
            cos_t = numpy.cos(theta)
            return numpy.array(
                [
                    +2 * a * (x[0] * cos_t + x[1] * sin_t) ** 2,
                    +2 * b * (x[0] * sin_t - x[1] * cos_t) ** 2,
                    +(a ** 2)
                    * 2
                    * (x[0] * cos_t + x[1] * sin_t)
                    * (-x[0] * sin_t + x[1] * cos_t)
                    + b ** 2
                    * 2
                    * (x[0] * sin_t - x[1] * cos_t)
                    * (x[0] * cos_t + x[1] * sin_t),
                ]
            ).T

        # We need to use some optimization here to find the new ellipses which best fit
        # the modified data. If xy_to_2d is the identity, we wouldn't need this.
        #
        # out = leastsq(f_ellipse, [1.0, 1.0, 0.0], Dfun=jac, full_output=True)
        # print(out)
        (a, b, theta), _ = leastsq(f_ellipse, [1.0, 1.0, 0.0], Dfun=jac)

        # (a, b, theta), _, infodict, msg, ierr = \
        #     leastsq(f, [1.0, 1.0, 0.0], full_output=True, Dfun=jac)
        # print(infodict['nfev'])

        # plot the scaled ellipse
        e = Ellipse(
            xy=tcenter,
            width=ellipse_scaling * 2 / a,
            height=ellipse_scaling * 2 / b,
            angle=theta / numpy.pi * 180,
            label=label,
        )
        ax.add_artist(e)
        e.set_alpha(alpha)
        e.set_facecolor(facecolor)
    return


def show_straights(cs):
    from mpl_toolkits.mplot3d import Axes3D

    # Some straight lines in XYZ
    t = numpy.linspace(0.0, 1.0, 101)
    n = 10

    fig = plt.figure()
    ax = fig.gca(projection=Axes3D.name)
    # ax.set_aspect('equal')

    for _ in range(n):
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


def xy_gamut_mesh(lcar):
    import optimesh
    import pygmsh

    observer = observers.cie_1931_2()

    # Gather all points on the horseshoe outline
    lmbda = 1.0e-9 * numpy.arange(380, 701)
    all_points = numpy.empty((len(lmbda), 2))
    for k in range(len(lmbda)):
        data = numpy.zeros(len(lmbda))
        data[k] = 1.0
        all_points[k] = _xyy_from_xyz100(spectrum_to_xyz100((lmbda, data), observer))[
            :2
        ]

    # Generate gmsh geometry: spline + straight line
    all_points = numpy.column_stack([all_points, numpy.zeros(len(all_points))])
    geom = pygmsh.built_in.Geometry()
    gmsh_points = [geom.add_point(pt, lcar) for pt in all_points]
    s1 = geom.add_spline(gmsh_points)
    s2 = geom.add_line(gmsh_points[-1], gmsh_points[0])
    ll = geom.add_line_loop([s1, s2])
    geom.add_plane_surface(ll)

    mesh = pygmsh.generate_mesh(geom)
    points, cells = optimesh.cvt.quasi_newton_uniform_lloyd(
        mesh.points, mesh.cells["triangle"], 1.0e-2, 100, omega=2.0
    )
    return points, cells


def get_mono_outline_xy(observer, max_stepsize, max_angle=None):
    """Monochromatic light of different frequencies form a horseshoe-like shape in
    xy-space. Get the outline of that space.
    """
    lmbda, data = observer

    m = lmbda.shape[0]
    mono = numpy.zeros(m)

    # first the straight connector at the bottom
    mono[:] = 0.0
    mono[-1] = 1.0
    first = _xyy_from_xyz100(spectrum_to_xyz100((lmbda, mono), observer))[:2]
    mono[:] = 0.0
    mono[0] = 1.0
    last = _xyy_from_xyz100(spectrum_to_xyz100((lmbda, mono), observer))[:2]
    #
    diff = first - last
    dist = numpy.sqrt(numpy.sum(diff ** 2))
    num_steps = dist / max_stepsize
    num_steps = int(num_steps) + 2
    # connection between lowest and highest frequencies
    vals_conn = numpy.array(
        [first * (1 - t) + last * t for t in numpy.linspace(0, 1, num_steps)]
    )

    vals_mono = [vals_conn[-1]]
    for k in range(1, m):
        mono[:] = 0.0
        mono[k] = 1.0
        val = _xyy_from_xyz100(spectrum_to_xyz100((lmbda, mono), observer))[:2]

        diff = vals_mono[-1] - val
        dist = numpy.sqrt(numpy.dot(diff, diff))

        if dist > max_stepsize:
            vals_mono.append(val)
    vals_mono.append(vals_conn[0])
    vals_mono = numpy.array(vals_mono)

    return vals_mono, vals_conn
