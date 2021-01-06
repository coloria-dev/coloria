import numpy
import pytest

import colorio


def test_flat_gamut():
    colorio.show_flat_gamut()


@pytest.mark.parametrize(
    "a", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion_variants(a):
    b = a + 1.0e-3 * numpy.random.rand(*a.shape)
    diff = colorio.delta(a, b)
    assert diff.shape == a.shape[1:]


def test_xy_gamut_mesh():
    points, cells = colorio.xy_gamut_mesh(0.05)

    # import meshio
    # meshio.write_points_cells("test.vtu", points, {"triangle": cells})
