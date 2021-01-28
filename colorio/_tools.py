import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from . import observers
from ._helpers import _find_Y
from .cs import HdrLinear, SrgbLinear
from .illuminants import planckian_radiator, spectrum_to_xyz100


def _xyy_from_xyz100(xyz):
    sum_xyz = np.sum(xyz, axis=0)
    x = xyz[0]
    y = xyz[1]
    return np.array([x / sum_xyz, y / sum_xyz, y / 100])


def _plot_monochromatic(observer, xy_to_2d, fill_horseshoe=True):
    # draw outline of monochromatic spectra
    lmbda = 1.0e-9 * np.arange(380, 701)
    values = []
    # TODO vectorize (see <https://github.com/numpy/numpy/issues/10439>)
    for k, _ in enumerate(lmbda):
        data = np.zeros(len(lmbda))
        data[k] = 1.0
        values.append(_xyy_from_xyz100(spectrum_to_xyz100((lmbda, data), observer))[:2])
    values = np.array(values)

    # Add the values between the first and the last point of the horseshoe
    t = np.linspace(0.0, 1.0, 101)
    connect = xy_to_2d(np.outer(values[0], t) + np.outer(values[-1], 1 - t))
    values = xy_to_2d(values.T).T
    full = np.concatenate([values, connect.T])

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
    for temp in np.arange(1000, 20001, 100):
        xyy_vals = xy_to_2d(
            _xyy_from_xyz100(spectrum_to_xyz100(planckian_radiator(temp), observer))
        )
        values.append(xyy_vals)
    values = np.array(values)
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
    lmbda = 1.0e-9 * np.arange(380, 701)
    all_points = np.empty((len(lmbda), 2))
    for k in range(len(lmbda)):
        data = np.zeros(len(lmbda))
        data[k] = 1.0
        all_points[k] = _xyy_from_xyz100(spectrum_to_xyz100((lmbda, data), observer))[
            :2
        ]

    # Generate gmsh geometry: spline + straight line
    all_points = np.column_stack([all_points, np.zeros(len(all_points))])
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
    mono = np.zeros(m)

    # first the straight connector at the bottom
    mono[:] = 0.0
    mono[-1] = 1.0
    first = _xyy_from_xyz100(spectrum_to_xyz100((lmbda, mono), observer))[:2]
    mono[:] = 0.0
    mono[0] = 1.0
    last = _xyy_from_xyz100(spectrum_to_xyz100((lmbda, mono), observer))[:2]
    #
    diff = first - last
    dist = np.sqrt(np.sum(diff ** 2))
    num_steps = dist / max_stepsize
    num_steps = int(num_steps) + 2
    # connection between lowest and highest frequencies
    vals_conn = np.array(
        [first * (1 - t) + last * t for t in np.linspace(0, 1, num_steps)]
    )

    vals_mono = [vals_conn[-1]]
    for k in range(1, m):
        mono[:] = 0.0
        mono[k] = 1.0
        val = _xyy_from_xyz100(spectrum_to_xyz100((lmbda, mono), observer))[:2]

        diff = vals_mono[-1] - val
        dist = np.sqrt(np.dot(diff, diff))

        if dist > max_stepsize:
            vals_mono.append(val)
    vals_mono.append(vals_conn[0])
    vals_mono = np.array(vals_mono)

    return vals_mono, vals_conn


def save_visible_gamut(colorspace, observer, illuminant, filename):
    import meshio
    from scipy.spatial import ConvexHull

    lmbda, illu = illuminant
    values = []

    # Iterate over every possible illuminant input and store it in values
    n = len(lmbda)
    values = np.empty((n * (n - 1) + 2, 3))
    k = 0

    # No light
    data = np.zeros(n)
    values[k] = spectrum_to_xyz100((lmbda, illu * data), observer=observer)
    k += 1
    # frequency blocks
    for width in range(1, n):
        data = np.zeros(n)
        data[:width] = 1.0
        for _ in range(n):
            values[k] = spectrum_to_xyz100((lmbda, illu * data), observer=observer)
            k += 1
            data = np.roll(data, shift=1)
    # Full illuminant
    data = np.ones(len(lmbda))
    values[k] = spectrum_to_xyz100((lmbda, illu * data), observer=observer)
    k += 1

    # scale the values such that the Y-coordinate of the white point (last entry)
    # has value 100.
    values *= 100 / values[-1][1]

    cells = ConvexHull(values).simplices

    if not colorspace.is_origin_well_defined:
        values = values[1:]
        cells = cells[~np.any(cells == 0, axis=1)]
        cells -= 1

    pts = colorspace.from_xyz100(values.T).T

    meshio.write_points_cells(filename, pts, cells={"triangle": cells})


def save_rgb_gamut(colorspace, filename: str, variant: str = "srgb", n: int = 50):
    import meshio
    import meshzoo

    if variant.lower() in ["srgb", "rec709"]:
        rgb_linear = SrgbLinear()
    else:
        assert variant.lower() in ["hdr", "rec2020", "rec2100"]
        rgb_linear = HdrLinear()

    points, cells = meshzoo.cube(nx=n, ny=n, nz=n)

    if not colorspace.is_origin_well_defined:
        # cut off [0, 0, 0] to avoid division by 0 in the xyz conversion
        points = points[1:]
        cells = cells[~np.any(cells == 0, axis=1)]
        cells -= 1

    pts = colorspace.from_rgb_linear(points.T).T
    # pts = colorspace.from_xyz100(rgb_linear.to_xyz100(points.T)).T
    assert pts.shape[1] == 3
    rgb = rgb_linear.to_rgb1(points)
    meshio.write_points_cells(filename, pts, {"tetra": cells}, point_data={"srgb": rgb})


def save_cone_gamut(colorspace, filename, observer, max_Y):
    import meshio
    import pygmsh

    with pygmsh.geo.Geometry() as geom:
        max_stepsize = 4.0e-2
        xy, _ = get_mono_outline_xy(observer, max_stepsize=max_stepsize)

        # append third component
        xy = np.column_stack([xy, np.full(xy.shape[0], 1.0e-5)])

        # Draw a cross.
        poly = geom.add_polygon(xy, mesh_size=max_stepsize)

        axis = [0, 0, max_Y]

        geom.extrude(poly, translation_axis=axis)

        mesh = geom.generate_mesh(verbose=False)
    # meshio.write(filename, mesh)

    pts = colorspace.from_xyz100(_xyy_to_xyz100(mesh.points.T)).T
    meshio.write_points_cells(filename, pts, {"tetra": mesh.get_cells_type("tetra")})


def show_visible_slice(*args, **kwargs):
    plt.figure()
    plot_visible_slice(*args, **kwargs)
    plt.show()
    plt.close()


def save_visible_slice(filename, *args, **kwargs):
    plt.figure()
    plot_visible_slice(*args, **kwargs)
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


def plot_visible_slice(colorspace, lightness, outline_prec=1.0e-2, fill_color="0.8"):
    # first plot the monochromatic outline
    mono_xy, conn_xy = get_mono_outline_xy(
        observer=observers.cie_1931_2(), max_stepsize=outline_prec
    )

    mono_vals = np.array([_find_Y(colorspace, xy, lightness) for xy in mono_xy])
    conn_vals = np.array([_find_Y(colorspace, xy, lightness) for xy in conn_xy])

    k1, k2 = [k for k in [0, 1, 2] if k != colorspace.k0]
    plt.plot(mono_vals[:, k1], mono_vals[:, k2], "-", color="k")
    plt.plot(conn_vals[:, k1], conn_vals[:, k2], ":", color="k")
    #
    if fill_color is not None:
        xyz = np.vstack([mono_vals, conn_vals[1:]])
        plt.fill(xyz[:, k1], xyz[:, k2], facecolor=fill_color, zorder=0)

    plt.axis("equal")
    plt.xlabel(colorspace.labels[k1])
    plt.ylabel(colorspace.labels[k2])
    plt.title(f"{colorspace.labels[colorspace.k0]} = {lightness}")


def show_rgb_slice(*args, **kwargs):
    plt.figure()
    plot_rgb_slice(*args, **kwargs)
    plt.show()
    plt.close()


def save_rgb_slice(filename, *args, **kwargs):
    plt.figure()
    plot_rgb_slice(*args, **kwargs)
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


def plot_rgb_slice(colorspace, lightness: float, n: int = 50, variant: str = "srgb"):
    import meshzoo

    # TODO HDR
    assert variant in ["srgb", "rec709"]

    # Get all RGB values that sum up to 1.
    srgb_vals, triangles = meshzoo.triangle(n=n)
    srgb_vals = srgb_vals.T

    # Use bisection to
    srgb_linear = SrgbLinear()
    tol = 1.0e-5
    # Use zeros() instead of empty() here to avoid invalid values when setting up
    # the cmap below.
    colorspace_vals = np.zeros((srgb_vals.shape[0], 3))
    srgb_linear_vals = np.zeros((srgb_vals.shape[0], 3))
    mask = np.ones(srgb_vals.shape[0], dtype=bool)
    for k, val in enumerate(srgb_vals):
        alpha_min = 0.0
        xyz100 = srgb_linear.to_xyz100(val * alpha_min)
        colorspace_val_min = colorspace.from_xyz100(xyz100)[colorspace.k0]
        if colorspace_val_min > lightness:
            mask[k] = False
            continue

        alpha_max = 1.0 / np.max(val)

        xyz100 = srgb_linear.to_xyz100(val * alpha_max)
        colorspace_val_max = colorspace.from_xyz100(xyz100)[colorspace.k0]
        if colorspace_val_max < lightness:
            mask[k] = False
            continue

        # bisection
        while True:
            alpha = (alpha_max + alpha_min) / 2
            srgb_linear_vals[k] = val * alpha
            xyz100 = srgb_linear.to_xyz100(srgb_linear_vals[k])
            colorspace_val = colorspace.from_xyz100(xyz100)
            if abs(colorspace_val[colorspace.k0] - lightness) < tol:
                break
            elif colorspace_val[colorspace.k0] < lightness:
                alpha_min = alpha
            else:
                assert colorspace_val[colorspace.k0] > lightness
                alpha_max = alpha
        colorspace_vals[k] = colorspace_val

    # Remove all triangles which have masked corner points
    tri_mask = np.all(mask[triangles], axis=1)
    if ~np.any(tri_mask):
        return
    triangles = triangles[tri_mask]

    # Unfortunately, one cannot use tripcolors with explicit RGB specification (see
    # <https://github.com/matplotlib/matplotlib/issues/10265>). As a workaround,
    # associate range(n) data with the points and create a colormap that associates
    # the integer values with the respective RGBs.
    z = np.arange(srgb_vals.shape[0])
    rgb = srgb_linear.to_rgb1(srgb_linear_vals)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("gamut", rgb, N=len(rgb))

    k1, k2 = [k for k in [0, 1, 2] if k != colorspace.k0]

    plt.tripcolor(
        colorspace_vals[:, k1],
        colorspace_vals[:, k2],
        triangles,
        z,
        shading="gouraud",
        cmap=cmap,
    )
    # plt.triplot(colorspace_vals[:, k1], colorspace_vals[:, k2], triangles=triangles)


def show_srgb1_gradient(*args, **kwargs):
    plt.figure()
    plot_srgb1_gradient(*args, **kwargs)
    plt.show()
    plt.close()


def plot_srgb1_gradient(colorspace, srgb0, srgb1, n=256):
    srgb = get_srgb1_gradient(colorspace, srgb0, srgb1, n=n)

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("empty", srgb, n)

    gradient = np.linspace(0.0, 1.0, n)
    gradient = np.vstack((gradient, gradient))
    plt.imshow(gradient, aspect="auto", cmap=cmap)
    plt.axis("off")
    plt.title(f"SRGB gradient in {colorspace.name}")


def get_srgb1_gradient(colorspace, srgb0, srgb1, n):
    # convert to colorspace
    cs = [colorspace.from_rgb1(srgb0), colorspace.from_rgb1(srgb1)]

    # linspace
    ls = np.linspace(cs[0], cs[1], endpoint=True, num=n, axis=0)

    # back to srgb
    srgb = colorspace.to_rgb1(ls.T).T

    srgb[srgb < 0] = 0.0
    srgb[srgb > 1] = 1.0
    return srgb


def show_srgb255_gradient(colorspace, srgb0, srgb1, n=256):
    srgb0 = np.asarray(srgb0)
    srgb1 = np.asarray(srgb1)
    show_srgb1_gradient(colorspace, srgb0 / 255, srgb1 / 255, n)


def plot_srgb255_gradient(colorspace, srgb0, srgb1, n=256):
    srgb0 = np.asarray(srgb0)
    srgb1 = np.asarray(srgb1)
    plot_srgb1_gradient(colorspace, srgb0 / 255, srgb1 / 255, n)


def get_srgb255_gradient(colorspace, srgb0, srgb1, n):
    srgb0 = np.asarray(srgb0)
    srgb1 = np.asarray(srgb1)
    return get_srgb1_gradient(colorspace, srgb0 / 255, srgb1 / 255, n) * 255


def show_primary_srgb_gradients(*args, **kwargs):
    plot_primary_srgb_gradients(*args, **kwargs)
    plt.show()
    plt.close()


def plot_primary_srgb_gradients(colorspace, n=256):
    pairs = [
        [([1, 1, 1], [1, 0, 0]), ([1, 0, 0], [0, 1, 0])],
        [([1, 1, 1], [0, 1, 0]), ([0, 1, 0], [0, 0, 1])],
        [([1, 1, 1], [0, 0, 1]), ([0, 0, 1], [1, 0, 0])],
        [([0, 0, 0], [1, 0, 0]), ([1, 0, 0], [0, 1, 1])],
        [([0, 0, 0], [0, 1, 0]), ([0, 1, 0], [1, 0, 1])],
        [([0, 0, 0], [0, 0, 1]), ([0, 0, 1], [1, 1, 0])],
    ]
    fig, axes = plt.subplots(len(pairs), 2)
    for i in range(len(pairs)):
        for j in range(2):
            pair = pairs[i][j]
            ax = axes[i][j]
            srgb = get_srgb1_gradient(colorspace, pair[0], pair[1], n=n)

            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", srgb, n)

            gradient = np.linspace(0.0, 1.0, n)
            gradient = np.vstack((gradient, gradient))
            ax.imshow(gradient, aspect="auto", cmap=cmap)
            ax.axis("off")
    fig.suptitle(f"primary SRGB gradients in {colorspace.name}")


def _xyy_to_xyz100(xyy):
    x, y, Y = xyy
    return np.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100
