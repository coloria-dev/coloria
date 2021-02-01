import matplotlib.pyplot as plt
import numpy as np

from . import observers
from ._helpers import _find_Y
from ._tools import get_mono_outline_xy, spectrum_to_xyz100
from .cs import XYY


def _get_visible_gamut_mesh(observer, max_Y1, h=4.0e-2):
    import meshio
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


def save_visible_gamut(filename, colorspace, observer, max_Y1, h=4.0e-2):
    import meshio

    points, cells = _get_visible_gamut_mesh(observer, max_Y1, h=h)

    xyz100 = XYY(1).to_xyz100(points.T)
    xyz100[xyz100 < 0] = 0.0
    points = colorspace.from_xyz100(xyz100).T

    meshio.write_points_cells(filename, points, {"tetra": cells})


def show_visible_gamut(colorspace, observer, max_Y1, show_grid=True, h=4.0e-2):
    import pyvista as pv
    import vtk

    points, cells = _get_visible_gamut_mesh(observer, max_Y1, h=h)

    xyz100 = XYY(1).to_xyz100(points.T)
    xyz100[xyz100 < 0] = 0.0
    points = colorspace.from_xyz100(xyz100).T

    cells = np.column_stack(
        [np.full(cells.shape[0], cells.shape[1], dtype=cells.dtype), cells]
    )

    # each cell is a VTK_HEXAHEDRON
    celltypes = np.full(len(cells), vtk.VTK_TETRA)

    grid = pv.UnstructuredGrid(cells.ravel(), celltypes, points)
    # grid.plot()

    p = pv.Plotter()
    p.add_mesh(grid)
    if show_grid:
        p.show_grid(
            xlabel=colorspace.labels[0],
            ylabel=colorspace.labels[1],
            zlabel=colorspace.labels[2],
        )
    p.show()


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
    plt.title(
        f"visible gamut slice in {colorspace.name} with "
        f"{colorspace.labels[colorspace.k0]}={lightness}"
    )
