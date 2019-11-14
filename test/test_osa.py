import numpy
import pytest

import colorio

numpy.random.seed(0)


@pytest.mark.parametrize(
    "xyz",
    [
        100 * numpy.random.rand(3),
        100 * numpy.random.rand(3, 7),
        # 100 * numpy.random.rand(3, 4, 5),
    ],
)
def test_conversion(xyz):
    osa = colorio.Osa()
    out = osa.to_xyz100(osa.from_xyz100(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-10)
    return