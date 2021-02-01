import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from ._nonlinear import regula_falsi
from .cs import HdrLinear, SrgbLinear


def save_rgb_gamut(filename: str, colorspace, variant: str = "srgb", n: int = 50):
    import meshio
    import meshzoo

    if variant.lower() in ["srgb", "rec709"]:
        rgb_linear = SrgbLinear()
    else:
        assert variant.lower() in ["hdr", "rec2020", "rec2100"]
        rgb_linear = HdrLinear()

    points, cells = meshzoo.cube_hexa((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), n)

    if not colorspace.is_origin_well_defined:
        # cut off [0, 0, 0] to avoid division by 0 in the xyz conversion
        points = points[1:]
        cells = cells[~np.any(cells == 0, axis=1)]
        cells -= 1

    pts = colorspace.from_rgb_linear(points.T).T
    # pts = colorspace.from_xyz100(rgb_linear.to_xyz100(points.T)).T
    assert pts.shape[1] == 3
    rgb = rgb_linear.to_rgb1(points)
    meshio.write_points_cells(
        filename, pts, {"hexahedron": cells}, point_data={"srgb": rgb}
    )


def show_rgb_gamut(
    colorspace, n: int = 51, show_grid: bool = True, camera_position=None
):
    import meshzoo
    import pyvista as pv
    import vtk

    points, cells = meshzoo.cube_hexa((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), n)
    cells = np.column_stack([np.full(cells.shape[0], cells.shape[1]), cells])

    srgb_linear = SrgbLinear()
    xyz100_coords = srgb_linear.to_xyz100(points.T)
    cs_coords = colorspace.from_xyz100(xyz100_coords).T

    # each cell is a VTK_HEXAHEDRON
    celltypes = np.full(len(cells), vtk.VTK_HEXAHEDRON, dtype=np.uint8)

    grid = pv.UnstructuredGrid(cells.ravel(), celltypes, cs_coords)
    # grid = grid.slice_orthogonal()
    # grid.slice_along_axis(n=7, axis="z")
    # single_slice = mesh.slice(normal=[0, 0, 1])

    p = pv.Plotter()
    p.add_mesh(
        grid,
        scalars=srgb_linear.to_rgb1(points.T).T,
        rgb=True,
        # show_edges=True,
    )
    if show_grid:
        p.show_grid(
            xlabel=colorspace.labels[0],
            ylabel=colorspace.labels[1],
            zlabel=colorspace.labels[2],
        )
    # camera_location = (0.5, 2.0, -2.0)
    # focus_point = (0.5, 0.0, 0.0)
    # viewup_vector = [0.0, 0.0, 0.0]
    # viewup_vector[colorspace.k0] = 1.0
    if camera_position is not None:
        p.camera_position = camera_position
    last_camera_position = p.show()
    return last_camera_position


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


def plot_rgb_slice(
    colorspace,
    lightness: float,
    n: int = 51,
    variant: str = "srgb",
    tol: float = 1.0e-5,
):
    # The best way to produce a slice is via a dedicated software, e.g., VTK. To avoid
    # the dependency, we're doing the following here: Take a slice out of the sRGB cube
    # and  project it into the the colorspace lightness level. Cut off all triangles
    # which don't fit. This cutting leads to jagged edges of the slice, but all other
    # methods have their imperfections, too.
    import meshzoo

    # TODO HDR
    assert variant in ["srgb", "rec709"]

    # Get all RGB values that sum up to 1.
    srgb_vals, triangles = meshzoo.triangle(n=n)
    srgb_vals = srgb_vals.T

    srgb_linear = SrgbLinear()
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

        def f(alpha):
            srgb_linear_vals[k] = val * alpha
            xyz100 = srgb_linear.to_xyz100(srgb_linear_vals[k])
            colorspace_val = colorspace.from_xyz100(xyz100)
            return colorspace_val[colorspace.k0] - lightness

        a, b = regula_falsi(f, alpha_min, alpha_max, tol)
        a = (a + b) / 2

        srgb_linear_vals[k] = val * a
        xyz100 = srgb_linear.to_xyz100(srgb_linear_vals[k])
        colorspace_val = colorspace.from_xyz100(xyz100)
        colorspace_vals[k] = colorspace_val

    # import meshplex
    # meshplex.MeshTri(colorspace_vals[:, 1:], triangles).show(show_coedges=False)

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

    import dufte

    plt.style.use(dufte.style)
    # default margins:
    plt.gca().margins(0.05)

    plt.tripcolor(
        colorspace_vals[:, k1],
        colorspace_vals[:, k2],
        triangles,
        z,
        shading="gouraud",
        cmap=cmap,
    )
    # plt.triplot(colorspace_vals[:, k1], colorspace_vals[:, k2], triangles=triangles)
    plt.title(
        f"sRGB gamut slice in {colorspace.name} with "
        f"{colorspace.labels[colorspace.k0]}={lightness}"
    )
    plt.xlabel(colorspace.labels[k1])
    plt.ylabel(colorspace.labels[k2])
    plt.gca().set_aspect("equal")
