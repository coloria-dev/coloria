import npx
import numpy as np
import pytest

import colorio

rng = np.random.default_rng(0)

wp1 = colorio.illuminants.whitepoints_cie1931["D65"]
wp2 = colorio.illuminants.whitepoints_cie1964["C"]

transforms = [
    (colorio.cat.bianco_schettini(wp1, wp2), 1.0e-13),
    (colorio.cat.bianco_schettini_pos(wp1, wp2), 1.0e-13),
    (colorio.cat.bradford(wp1, wp2), 1.0e-13),
    (colorio.cat.cat02(wp1, wp2, 1.0, 20.0), 1.0e-13),
    (colorio.cat.cat16(wp1, wp2, 1.0, 20.0), 1.0e-13),
    (colorio.cat.cmccat2000(wp1, wp2, 1.0, 20.0, 30.0), 1.0e-13),
    (colorio.cat.sharp(wp1, wp2), 1.0e-13),
    (colorio.cat.von_kries(wp1, wp2), 1.0e-13),
]


@pytest.mark.parametrize("cat, tol", transforms)
def test_zero(cat, tol):
    print(cat)
    xyz = np.zeros(3)
    out = cat[1] @ (cat[0] @ xyz)
    print(xyz)
    print(out)
    assert np.all(np.abs(xyz - out) < tol)


@pytest.mark.parametrize("cat, tol", transforms)
@pytest.mark.parametrize(
    "xyz",
    [
        100 * rng.random(3),
        100 * rng.random((3, 7)),
        100 * rng.random((3, 4, 5)),
        # make sure the thing works with lists as input, too
        [1.0, 2.0, 3.0],
    ],
)
def test_round_trip(cat, tol, xyz):
    xyz_orig = xyz.copy()
    print(xyz)
    out = npx.dot(cat[1], npx.dot(cat[0], xyz))
    # make sure that xyz doesn't change during the calls
    assert np.all(np.abs(np.asarray(xyz) - xyz_orig) < tol * np.abs(xyz_orig))
    print(xyz)
    print(cat)
    print(out)
    print(np.max(np.abs(xyz - out)))
    assert np.asarray(xyz).shape == out.shape
    assert np.all(np.abs(np.asarray(xyz) - out) < tol * np.abs(xyz))

    # # make sure it works the other way around, too
    # out = cat.apply(cat.apply_inv(xyz))
    # assert np.all(np.abs(np.array(xyz) - out) < tol * np.abs(xyz))
