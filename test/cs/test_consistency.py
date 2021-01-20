import numpy as np
import pytest

import colorio

np.random.seed(0)

colorspaces = [
    (colorio.cs.CIELCH(), 1.0e-14),
    (colorio.cs.CIEHCL(), 1.0e-13),
    (colorio.cs.CIELAB(), 1.0e-14),
    (colorio.cs.CIELUV(), 1.0e-14),
    (colorio.cs.DIN99(), 1.0e-14),
    (colorio.cs.HdrLinear(), 1.0e-14),
    (colorio.cs.IPT(), 1.0e-14),
    (colorio.cs.JzAzBz(), 1.0e-12),
    (colorio.cs.OKLAB(), 1.0e-14),
    (colorio.cs.OsaUcs(), 1.0e-11),
    (colorio.cs.PROLAB(), 1.0e-14),
    (colorio.cs.RLAB(), 1.0e-14),
    (colorio.cs.XYY(1), 1.0e-14),
    (colorio.cs.XYY(100), 1.0e-14),
    (colorio.cs.XYZ(1), 1.0e-14),
    (colorio.cs.XYZ(100), 1.0e-14),
]


@pytest.mark.parametrize(
    "cs, tol", [(cs, tol) for cs, tol in colorspaces if cs.is_origin_well_defined]
)
def test_zero(cs, tol):
    print(cs)
    xyz = np.zeros(3)
    out = cs.to_xyz100(cs.from_xyz100(xyz))
    print(xyz)
    print(out)
    assert np.all(np.abs(xyz - out) < tol)


@pytest.mark.parametrize("cs, tol", colorspaces)
@pytest.mark.parametrize(
    "xyz",
    [
        100 * np.random.rand(3),
        100 * np.random.rand(3, 7),
        100 * np.random.rand(3, 4, 5),
        # make sure the thing works with lists and input, too
        [1.0, 2.0, 3.0],
    ],
)
def test_conversion(cs, tol, xyz):
    xyz_orig = xyz.copy()
    print(xyz)
    out = cs.to_xyz100(cs.from_xyz100(xyz))
    # make sure that xyz doesn't change during the calls
    assert np.all(np.abs(np.asarray(xyz) - xyz_orig) < tol * np.abs(xyz_orig))
    print(xyz)
    print(cs)
    print(out)
    assert np.all(np.abs(np.array(xyz) - out) < tol * np.abs(xyz))

    # # make sure it works the other way around, too
    # out = cs.from_xyz100(cs.to_xyz100(xyz))
    # assert np.all(np.abs(np.array(xyz) - out) < tol * np.abs(xyz))
