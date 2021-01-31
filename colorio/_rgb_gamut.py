import matplotlib.pyplot as plt
import numpy as np

from .cs import HdrLinear, SrgbLinear


def save_rgb_gamut(colorspace, filename: str, variant: str = "srgb", n: int = 50):
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


def show_rgb_gamut(colorspace, n: int = 51):
    import meshzoo
    import pyvista as pv
    import vtk

    points, cells = meshzoo.cube_hexa((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), n)
    cells = np.column_stack([np.full(cells.shape[0], 8), cells])

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
    mesh = p.add_mesh(
        grid,
        scalars=srgb_linear.to_rgb1(points.T).T,
        rgb=True,
        # show_edges=True,
    )
    p.show()


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
    import meshzoo
    from scipy.spatial import Delaunay

    # TODO HDR
    assert variant in ["srgb", "rec709"]

    srgb_linear = SrgbLinear()

    def f(pts):
        # transform to colorspace coordinates
        xyz100_coords = srgb_linear.to_xyz100(pts)
        cs_coords = colorspace.from_xyz100(xyz100_coords)
        return cs_coords[colorspace.k0] - lightness

    srgb_points, cells = meshzoo.cube_hexa((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), n)
    is_lightness_above = f(srgb_points.T) > 0

    cells_ila = is_lightness_above[cells]

    # find all cells where at least one corner point is below and at least one is above
    # the target lightness
    is_intermediate_cell = np.any(cells_ila, axis=1) & np.any(~cells_ila, axis=1)

    cells = cells[is_intermediate_cell]
    cells_ila = cells_ila[is_intermediate_cell]

    # Now collect all edges along which the target lightness is crossed
    local_edges = [
        (0, 1),
        (0, 3),
        (0, 4),
        (1, 2),
        (1, 5),
        (2, 3),
        (2, 6),
        (3, 7),
        (4, 5),
        (4, 7),
        (5, 6),
        (6, 7),
    ]
    edges = []
    for e0, e1 in local_edges:
        a = np.logical_xor(cells_ila[:, e0], cells_ila[:, e1])
        edges.append(np.column_stack([cells[a, e0], cells[a, e1]]))

    edges = np.concatenate(edges)
    print(edges)
    print()

    # make unique
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges, axis=0)

    # now bisect to identify the points with the target lightness
    print()
    print(srgb_points.shape)
    print()
    edge_pts = srgb_points[edges]

    a = edge_pts[:, 0]
    b = edge_pts[:, 1]

    print(a)
    print()
    print(b)
    print()

    # a, b = bisect(f, a.T, b.T, tol=1.0e-5)
    a, b = regula_falsi(f, a.T, b.T, tol=1.0e-5)
    sol = (a + b) / 2
    print(sol.T)
    print(lightness, colorspace.k0)
    print()

    xyz100_coords = srgb_linear.to_xyz100(sol)
    cs_coords = colorspace.from_xyz100(xyz100_coords)
    print(cs_coords.T)
    print(cs_coords.shape)

    # The cells could actually be derived from the connectivity info of the transformed
    # cube. Since this is getting a little messy, rather get a triangulation of the
    # convex hull and kick out all cells the barycenters of which aren't in the set.
    cells = Delaunay(cs_coords[1:].T).simplices
    barycenters = np.sum(cs_coords[1:, cells], axis=-1).T / 3
    barycenters = np.column_stack(
        [np.full(barycenters.shape[0], lightness), barycenters]
    )
    xyz100 = colorspace.to_xyz100(barycenters.T)
    srgb_vals = srgb_linear.from_xyz100(xyz100)
    tol = 1.0e-5
    is_legal_srgb = np.all(srgb_vals < 1.0 + tol, axis=0) & np.all(
        srgb_vals > -tol, axis=0
    )

    # print(np.where(~is_legal_srgb))
    # print(srgb_vals[:, 0])
    # print(cells[0])
    # print(cs_coords[:, 1273])
    # print(cs_coords[:, 75])
    # print(cs_coords[:, 37])
    # print()
    # print(srgb_linear.from_xyz100(colorspace.to_xyz100(cs_coords[:, 1273])))
    # print(srgb_linear.from_xyz100(colorspace.to_xyz100(cs_coords[:, 75])))
    # print(srgb_linear.from_xyz100(colorspace.to_xyz100(cs_coords[:, 37])))
    # exit(1)

    # cells = cells[is_legal_srgb]

    print(is_legal_srgb.shape, is_legal_srgb.sum())

    for cell in cells:
        for i0, i1 in [[0, 1], [1, 2], [2, 1]]:
            p0 = cs_coords[1:, cell[i0]]
            p1 = cs_coords[1:, cell[i1]]
            plt.plot([p0[0], p1[0]], [p0[1], p1[1]], "k-")

    for x in barycenters[~is_legal_srgb]:
        plt.plot(x[1], x[2], "x", color="C0")

    plt.scatter(*cs_coords[1:], marker="x", color="k")
    plt.gca().set_aspect("equal")
    plt.show()
    exit(1)

    print(edges)

    print(is_lightness_above[edges])
    exit(1)

    print(len(edges))
    exit(1)
    print(cells.shape)
    print(cells_ila.shape)

    exit(1)

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
