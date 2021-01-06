import numpy
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz",
    [
        numpy.random.rand(3) * 100,
        numpy.random.rand(3, 7) * 100,
        numpy.random.rand(3, 4, 5) * 100,
    ],
)
def test_conversion(xyz):
    ciehcl = colorio.cs.CIEHCL()
    out = ciehcl.to_xyz100(ciehcl.from_xyz100(xyz))
    assert numpy.all(numpy.abs(xyz - out) < 1.0e-13 * numpy.abs(xyz))
