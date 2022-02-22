from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike

from . import observers
from ._helpers import SpectralData
from .cs import ColorCoordinates, ColorSpace, convert, string_to_cs
from .illuminants import planckian_radiator, spectrum_to_xyz100


def _unit_vector(n, k):
    vec = np.zeros(n)
    vec[k] = 1.0
    return vec


def _plot_monochromatic(observer, fill_horseshoe=True):
    # draw outline of monochromatic spectra
    lmbda_nm = np.arange(380, 701)

    n = len(lmbda_nm)
    values = ColorCoordinates(
        np.array(
            [
                spectrum_to_xyz100(SpectralData(lmbda_nm, _unit_vector(n, k)), observer)
                for k in range(n)
            ]
        ).T,
        "XYZ100",
    )
    values = convert(values, "XYY100").hue.T

    # Add the values between the first and the last point of the horseshoe
    t = np.linspace(0.0, 1.0, 101)
    connect = np.outer(values[0], t) + np.outer(values[-1], 1 - t)
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


def _plot_planckian_locus(observer):
    # plot planckian locus
    values = ColorCoordinates(
        np.array(
            [
                spectrum_to_xyz100(planckian_radiator(temp), observer)
                for temp in np.arange(1000, 20001, 100)
            ]
        ).T,
        "XYZ100",
    )
    x, y = convert(values, "XYY100").data[:2]
    plt.plot(x, y, ":k", label="Planckian locus")


def plot_xy_gamut(fill_horseshoe=True, plot_planckian_locus=True):
    """Show a flat color gamut, by default xy. There exists a chroma gamut for all
    color models which transform lines in XYZ to lines, and hence have a natural
    decomposition into lightness and chroma components. Also, the flat gamut is the
    same for every lightness value. Examples for color models with this property are
    CIELUV and IPT, examples for color models without are CIELAB and CIECAM02.
    """
    observer = observers.cie_1931_2()
    # observer = observers.cie_1964_10()

    _plot_monochromatic(observer, fill_horseshoe=fill_horseshoe)
    # plt.grid()

    # if plot_rgb_triangle:
    #     _plot_rgb_triangle()
    if plot_planckian_locus:
        _plot_planckian_locus(observer)

    plt.gca().set_aspect("equal")
    # plt.legend()
    plt.xlabel("x")
    plt.ylabel("y")
    return plt


def xy_gamut_mesh(lcar):
    import optimesh
    import pygmsh

    observer = observers.cie_1931_2()

    # Gather all points on the horseshoe outline
    lmbda_nm = np.arange(380, 701)

    n = len(lmbda_nm)
    xyz100 = ColorCoordinates(
        np.array(
            [
                spectrum_to_xyz100(SpectralData(lmbda_nm, _unit_vector(n, k)), observer)
                for k in range(n)
            ]
        ).T,
        "XYZ100",
    )

    x, y = convert(xyz100, "XYY100").data[:2]

    # Generate gmsh geometry: spline + straight line
    all_points = np.column_stack([x, y, np.zeros_like(x)])
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
    m = observer.lmbda_nm.shape[0]
    mono = np.zeros(m)

    # first the straight connector at the bottom
    mono[:] = 0.0
    mono[-1] = 1.0
    mono_spectrum = SpectralData(observer.lmbda_nm, mono)
    first = convert(
        ColorCoordinates(spectrum_to_xyz100(mono_spectrum, observer), "XYZ100"),
        "XYY100",
    ).data[:2]

    mono[:] = 0.0
    mono[0] = 1.0
    mono_spectrum = SpectralData(observer.lmbda_nm, mono)
    last = convert(
        ColorCoordinates(spectrum_to_xyz100(mono_spectrum, observer), "XYZ100"),
        "XYY100",
    ).data[:2]

    #
    diff = first - last
    dist = np.sqrt(np.sum(diff**2))
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
        mono_spectrum = SpectralData(observer.lmbda_nm, mono)
        val = convert(
            ColorCoordinates(spectrum_to_xyz100(mono_spectrum, observer), "XYZ100"),
            "XYY100",
        ).data[:2]

        diff = vals_mono[-1] - val
        dist = np.sqrt(np.dot(diff, diff))

        if dist > max_stepsize:
            vals_mono.append(val)
    vals_mono.append(vals_conn[0])
    vals_mono = np.array(vals_mono)

    return vals_mono, vals_conn


def plot_srgb1_gradient(
    colorspace: ColorSpace | str, srgb0: ArrayLike, srgb1: ArrayLike, n: int = 256
):
    srgb = get_srgb1_gradient(colorspace, srgb0, srgb1, n=n)

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("empty", srgb, n)

    gradient = np.linspace(0.0, 1.0, n)
    gradient = np.vstack((gradient, gradient))
    plt.imshow(gradient, aspect="auto", cmap=cmap)
    plt.axis("off")
    plt.title(f"SRGB gradient in {colorspace.name}")
    return plt


def get_srgb1_gradient(
    colorspace: ColorSpace | str, srgb0: ArrayLike, srgb1: ArrayLike, n: int
) -> np.ndarray:
    # convert to colorspace
    cs0 = convert(ColorCoordinates(srgb0, "srgb_linear"), colorspace).data
    cs1 = convert(ColorCoordinates(srgb1, "srgb_linear"), colorspace).data

    # linspace
    ls = np.linspace(cs0, cs1, endpoint=True, num=n, axis=0)
    coords = ColorCoordinates(ls.T, colorspace)
    return convert(coords, "srgb1", mode="clip").data.T


def plot_srgb255_gradient(
    colorspace: ColorSpace | str, srgb0: ArrayLike, srgb1: ArrayLike, n: int = 256
):
    srgb0 = np.asarray(srgb0)
    srgb1 = np.asarray(srgb1)
    return plot_srgb1_gradient(colorspace, srgb0 / 255, srgb1 / 255, n)


def get_srgb255_gradient(
    colorspace: ColorSpace | str, srgb0: ArrayLike, srgb1: ArrayLike, n: int
) -> np.ndarray:
    srgb0 = np.asarray(srgb0)
    srgb1 = np.asarray(srgb1)
    return get_srgb1_gradient(colorspace, srgb0 / 255, srgb1 / 255, n) * 255


def plot_primary_srgb_gradients(colorspace: ColorSpace | str, n: int = 256):
    if isinstance(colorspace, str):
        colorspace = string_to_cs(colorspace)

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
    return plt
