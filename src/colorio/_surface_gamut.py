from __future__ import annotations

import numpy as np

from ._helpers import SpectralData
from ._tools import spectrum_to_xyz100
from .cs import ColorSpace, string_to_cs


def _get_surface_gamut_mesh(
    colorspace: ColorSpace | str, observer: SpectralData, illuminant: SpectralData
) -> tuple[np.ndarray, np.ndarray]:
    from scipy.spatial import ConvexHull

    if isinstance(colorspace, str):
        colorspace = string_to_cs(colorspace)

    # lmbda, illu = illuminant
    values = []

    # Iterate over every possible illuminant input and store it in values
    n = len(illuminant.lmbda_nm)
    mask = np.zeros(n)
    # frequency blocks
    values = []
    # no light
    xyz100 = spectrum_to_xyz100(SpectralData(illuminant.lmbda_nm, mask), observer)
    values.append(xyz100)
    for width in range(1, n):
        mask[:] = 0.0
        mask[:width] = 1.0
        for _ in range(n):
            spec = SpectralData(illuminant.lmbda_nm, mask * illuminant.data)
            xyz100 = spectrum_to_xyz100(spec, observer)
            assert np.all(xyz100 >= 0), xyz100
            values.append(xyz100)
            mask = np.roll(mask, shift=1)
    # Full illuminant
    values.append(spectrum_to_xyz100(illuminant, observer))
    values = np.array(values)

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


def plot_surface_gamut(
    colorspace: ColorSpace | str, observer, illuminant, show_grid=True
):
    import pyvista as pv
    import vtk

    if isinstance(colorspace, str):
        colorspace = string_to_cs(colorspace)

    points, cells = _get_surface_gamut_mesh(colorspace, observer, illuminant)
    cells = np.column_stack(
        [np.full(cells.shape[0], cells.shape[1], dtype=cells.dtype), cells]
    )

    # each cell is a VTK_HEXAHEDRON
    celltypes = np.full(len(cells), vtk.VTK_TRIANGLE)

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
    return p
