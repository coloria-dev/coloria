import matplotlib.pyplot as plt
import numpy as np

from . import observers
from ._helpers import _find_Y
from ._tools import get_mono_outline_xy, spectrum_to_xyz100


def _xyy_to_xyz100(xyy):
    x, y, Y = xyy
    return np.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100


def save_visible_gamut(colorspace, filename, observer, max_Y):
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
    plt.title(
        f"visible gamut slice in {colorspace.name} with "
        f"{colorspace.labels[colorspace.k0]}={lightness}"
    )
