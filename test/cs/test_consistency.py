import numpy
import pytest

import colorio

numpy.random.seed(0)


@pytest.mark.parametrize(
    "cs, tol",
    [
        (colorio.cs.CIELCH(), 1.0e-14),
        (colorio.cs.CIEHCL(), 1.0e-13),
        (colorio.cs.CIELAB(), 1.0e-14),
        (colorio.cs.CIELUV(), 1.0e-14),
        (colorio.cs.HdrLinear(), 1.0e-14),
        (colorio.cs.IPT(), 1.0e-14),
        (colorio.cs.JzAzBz(), 1.0e-12),
        (colorio.cs.OKLAB(), 1.0e-14),
        (colorio.cs.OsaUcs(), 1.0e-11),
        (colorio.cs.PROLAB(), 1.0e-14),
        (colorio.cs.RLAB(), 1.0e-14),
        (colorio.cs.XYY(1), 1.0e-14),
        (colorio.cs.XYY(100), 1.0e-14),
        (colorio.cs.XYZ1(), 1.0e-14),
    ],
)
@pytest.mark.parametrize(
    "xyz",
    [
        100 * numpy.random.rand(3),
        100 * numpy.random.rand(3, 7),
        100 * numpy.random.rand(3, 4, 5),
    ],
)
def test_conversion(cs, tol, xyz):
    xyz_orig = xyz.copy()
    print(xyz)
    out = cs.to_xyz100(cs.from_xyz100(xyz))
    # make sure that xyz doesn't change during the calls
    assert numpy.all(numpy.abs(xyz - xyz_orig) < tol * numpy.abs(xyz_orig))
    print(xyz)
    print(cs)
    print(out)
    assert numpy.all(numpy.abs(xyz - out) < tol * numpy.abs(xyz))
