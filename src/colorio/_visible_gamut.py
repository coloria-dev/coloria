from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from . import observers
from ._tools import get_mono_outline_xy
from .cs import ColorCoordinates, ColorSpace, convert, string_to_cs


def _get_visible_gamut_mesh(observer, max_Y1, h=4.0e-2):
    import pygmsh

    with pygmsh.geo.Geometry() as geom:
        xy, _ = get_mono_outline_xy(observer, max_stepsize=h)
        # append third component
        xy = np.column_stack([xy, np.full(xy.shape[0], 1.0e-5)])
        poly = geom.add_polygon(xy, mesh_size=h)
        axis = [0, 0, max_Y1]
        geom.extrude(poly, translation_axis=axis)
        mesh = geom.generate_mesh(verbose=False)

    return mesh.points, mesh.get_cells_type("tetra")


def plot_visible_gamut(
    colorspace: ColorSpace | str, observer, max_Y1, show_grid=True, h=4.0e-2
):
    import pyvista as pv
    import vtk

    if isinstance(colorspace, str):
        colorspace = string_to_cs(colorspace)

    points, cells = _get_visible_gamut_mesh(observer, max_Y1, h=h)

    xyy1 = ColorCoordinates(points.T, "XYY1")
    xyz100 = convert(xyy1, "XYZ100")
    xyz100.data[xyz100.data < 0] = 0.0
    coords = convert(xyz100, colorspace)

    # xyz100 = XYY(1).to_xyz100(points.T)
    # xyz100[xyz100 < 0] = 0.0
    # points = colorspace.from_xyz100(xyz100).T

    cells = np.column_stack(
        [np.full(cells.shape[0], cells.shape[1], dtype=cells.dtype), cells]
    )

    # each cell is a VTK_HEXAHEDRON
    celltypes = np.full(len(cells), vtk.VTK_TETRA)

    grid = pv.UnstructuredGrid(cells.ravel(), celltypes, coords.data.T)
    # grid.plot()

    p = pv.Plotter()
    p.add_mesh(grid)
    if show_grid:
        p.show_grid(
            xlabel=colorspace.labels[0],
            ylabel=colorspace.labels[1],
            zlabel=colorspace.labels[2],
        )

    return p


def plot_visible_slice(
    colorspace: ColorSpace | str,
    lightness: float,
    outline_prec: float = 1.0e-2,
    fill_color="0.8",
):
    if isinstance(colorspace, str):
        colorspace = string_to_cs(colorspace)

    # first plot the monochromatic outline
    mono_xy, conn_xy = get_mono_outline_xy(
        observer=observers.cie_1931_2(), max_stepsize=outline_prec
    )

    mono_vals = ColorCoordinates(
        np.array([_find_Y(colorspace, xy, lightness).data for xy in mono_xy]).T,
        colorspace,
    )
    conn_vals = ColorCoordinates(
        np.array([_find_Y(colorspace, xy, lightness).data for xy in conn_xy]).T,
        colorspace,
    )

    plt.plot(*mono_vals.hue, "-", color="k")
    plt.plot(*conn_vals.hue, ":", color="k")
    #
    if fill_color is not None:
        x, y = np.column_stack([mono_vals.hue, conn_vals.hue[:, 1:]])
        plt.fill(x, y, facecolor=fill_color, zorder=0)

    plt.axis("equal")
    plt.xlabel(colorspace.hue_labels[0])
    plt.ylabel(colorspace.hue_labels[1])
    plt.title(
        f"visible gamut slice in {colorspace.name} with "
        f"{colorspace.lightness_label}={lightness}"
    )

    return plt


def _find_Y(cs, xy, level, tol=1.0e-5):
    """Use bisection to find a matching Y value that projects the xy into the given
    level.
    """
    x, y = xy
    min_Y = 0.0

    xyy = ColorCoordinates([x, y, min_Y], "xyy1")
    min_val = convert(xyy, cs).lightness
    assert min_val <= level

    # search for an appropriate max_Y to start with
    max_Y = 1.0
    while convert(ColorCoordinates([x, y, max_Y], "xyy1"), cs).lightness < level:
        max_Y *= 2

    while True:
        Y = (max_Y + min_Y) / 2
        val = convert(ColorCoordinates([x, y, Y], "xyy1"), cs)
        if abs(val.lightness - level) < tol:
            break
        elif val.lightness > level:
            max_Y = Y
        else:
            assert val.lightness < level
            min_Y = Y

    return val
