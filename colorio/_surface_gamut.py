import numpy as np

from ._tools import spectrum_to_xyz100


def _get_surface_gamut_mesh(colorspace, observer, illuminant):
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
    return pts, cells


def save_surface_gamut(colorspace, observer, illuminant, filename):
    import meshio

    pts, cells = _get_surface_gamut_mesh(colorspace, observer, illuminant)
    meshio.write_points_cells(filename, pts, cells={"triangle": cells})


def show_surface_gamut(colorspace, observer, illuminant, show_grid=True):
    import pyvista as pv
    import vtk

    points, cells = _get_surface_gamut_mesh(colorspace, observer, illuminant)
    cells = np.column_stack([np.full(cells.shape[0], cells.shape[1]), cells])

    # each cell is a VTK_HEXAHEDRON
    celltypes = np.full(len(cells), vtk.VTK_TRIANGLE)

    grid = pv.UnstructuredGrid(cells.ravel(), celltypes, points, xlabel="Easting")
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
