import numpy as np
import pytest

import colorio

rng = np.random.default_rng(0)


@pytest.mark.parametrize(
    "xyz",
    [
        100 * rng.random(3),
        100 * rng.random((3, 7)),
        100 * rng.random((3, 4, 5)),
    ],
)
def test_conversion(xyz, tol=1.0e-12):
    # test with srgb conditions
    cam16 = colorio.cs.CAM16(0.69, 20, 20)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    out = cam16.to_xyz100(np.array([J, C, H]), "JCH")
    assert np.all(np.abs(xyz - out) < tol * np.abs(xyz))

    out = cam16.to_xyz100(np.array([Q, M, h]), "QMh")
    assert np.all(np.abs(xyz - out) < tol * np.abs(xyz))

    out = cam16.to_xyz100(np.array([J, s, h]), "Jsh")
    assert np.all(np.abs(xyz - out) < tol * np.abs(xyz))


@pytest.mark.parametrize("xyz", [np.zeros(3), np.zeros((3, 4, 5))])
def test_zero(xyz, tol=1.0e-12):
    cam16 = colorio.cs.CAM16(0.69, 20, 20)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    assert np.all(J == 0.0)
    assert np.all(C == 0.0)
    assert np.all(h == 0.0)
    assert np.all(M == 0.0)
    assert np.all(s == 0.0)
    assert np.all(Q == 0.0)

    out = cam16.to_xyz100(np.array([J, C, H]), "JCH")
    assert np.all(np.abs(out) < tol)

    out = cam16.to_xyz100(np.array([Q, M, h]), "QMh")
    assert np.all(np.abs(out) < tol)

    out = cam16.to_xyz100(np.array([J, s, h]), "Jsh")
    assert np.all(np.abs(out) < tol)


@pytest.mark.parametrize(
    "xyz,ref",
    [
        (
            [1.0, 0.0, 0.0],
            [2.2994065595734066, 113.32448472150614, 2.711228540807689],
        )
    ],
)
def test_reference_values(xyz, ref):
    cam16 = colorio.cs.CAM16UCS(0.69, 20, 20)
    out = cam16.from_xyz100(xyz)
    ref = np.array(ref)
    print(list(out))
    assert np.all(np.abs(ref - out) < 1.0e-14 * ref)


def test_whitepoint():
    # With infinite luminance of the adapting field, the whitepoint is found
    # at (100, 0, 0).
    L_A = np.inf
    cam16 = colorio.cs.CAM16UCS(0.69, 20, L_A)
    out = cam16.from_xyz100(colorio.illuminants.whitepoints_cie1931["D65"])
    assert np.all(out == [100, 0, 0])


if __name__ == "__main__":
    test_conversion(np.random.rand(3))
