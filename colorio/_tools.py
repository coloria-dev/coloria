import matplotlib.pyplot as plt
import numpy

from . import observers
from .illuminants import planckian_radiator, spectrum_to_xyz100


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


def show_flat_gamut(*args, **kwargs):
    plot_flat_gamut(*args, **kwargs)
    plt.show()


def plot_flat_gamut(
    xy_to_2d=lambda xy: xy,
    axes_labels=("x", "y"),
    # plot_rgb_triangle=True,
    fill_horseshoe=True,
    plot_planckian_locus=True,
):
    """Show a flat color gamut, by default xy. There exists a chroma gamut for all
    color models which transform lines in XYZ to lines, and hence have a natural
    decomposition into lightness and chroma components. Also, the flat gamut is the
    same for every lightness value. Examples for color models with this property are
    CIELUV and IPT, examples for color models without are CIELAB and CIECAM02.
    """
    observer = observers.cie_1931_2()
    # observer = observers.cie_1964_10()

    _plot_monochromatic(observer, xy_to_2d, fill_horseshoe=fill_horseshoe)
    # plt.grid()

    # if plot_rgb_triangle:
    #     _plot_rgb_triangle(xy_to_2d)
    if plot_planckian_locus:
        _plot_planckian_locus(observer, xy_to_2d)

    plt.gca().set_aspect("equal")
    # plt.legend()
    plt.xlabel(axes_labels[0])
    plt.ylabel(axes_labels[1])


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
    with pygmsh.geo.Geometry() as geom:
        gmsh_points = [geom.add_point(pt, lcar) for pt in all_points]
        s1 = geom.add_spline(gmsh_points)
        s2 = geom.add_line(gmsh_points[-1], gmsh_points[0])
        ll = geom.add_curve_loop([s1, s2])
        geom.add_plane_surface(ll)
        mesh = geom.generate_mesh()

    # Work around numpy bug <https://github.com/numpy/numpy/issues/17760>
    cells = mesh.get_cells_type("triangle").astype(int)
    points, cells = optimesh.optimize_points_cells(
        mesh.points, cells, "lloyd", 1.0e-2, 100, omega=2.0
    )
    return points, cells


def get_mono_outline_xy(observer, max_stepsize):
    """Monochromatic light of different frequencies form a horseshoe-like shape in
    xy-space. Get the outline of that space.
    """
    lmbda, _ = observer

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
