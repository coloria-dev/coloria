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
def test_conversion(xyz):
    # test with srgb conditions
    L_A = 64 / np.pi / 5
    cam16 = colorio.cs.CAM16(0.69, 20, L_A)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    out = cam16.to_xyz100(np.array([J, C, H]), "JCH")
    assert np.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = cam16.to_xyz100(np.array([Q, M, h]), "QMh")
    assert np.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = cam16.to_xyz100(np.array([J, s, h]), "Jsh")
    assert np.all(abs(xyz - out) < 1.0e-13 * abs(xyz))


@pytest.mark.parametrize("xyz", [np.zeros(3), np.zeros((3, 4, 5))])
def test_zero(xyz):
    cam16 = colorio.cs.CAM16(0.69, 20, 64 / np.pi / 5)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    assert np.all(J == 0.0)
    assert np.all(C == 0.0)
    assert np.all(h == 0.0)
    assert np.all(M == 0.0)
    assert np.all(s == 0.0)
    assert np.all(Q == 0.0)

    out = cam16.to_xyz100(np.array([J, C, H]), "JCH")
    assert np.all(abs(out) < 1.0e-13)

    out = cam16.to_xyz100(np.array([Q, M, h]), "QMh")
    assert np.all(abs(out) < 1.0e-13)

    out = cam16.to_xyz100(np.array([J, s, h]), "Jsh")
    assert np.all(abs(out) < 1.0e-13)


@pytest.mark.parametrize(
    "xyz", [np.random.rand(3), np.random.rand(3, 7), np.random.rand(3, 4, 5)]
)
def test_conversion_variants(xyz):
    # test with srgb conditions
    L_A = 64 / np.pi / 5
    cam16 = colorio.cs.CAM16UCS(0.69, 20, L_A)
    out = cam16.to_xyz100(cam16.from_xyz100(xyz))
    assert np.all(abs(xyz - out) < 1.0e-14)


@pytest.mark.parametrize(
    "xyz,ref",
    [
        (
            [1.0, 0.0, 0.0],
            [2.28402560268459e00, 1.01502029350636e02, 2.42718425228025e00],
        )
    ],
)
def test_reference_values(xyz, ref):
    L_A = 64 / np.pi / 5
    cam16 = colorio.cs.CAM16UCS(0.69, 20, L_A)
    out = cam16.from_xyz100(xyz)
    ref = np.array(ref)
    assert np.all(abs(ref - out) < 1.0e-14 * ref)


def test_whitepoint():
    # With infinite luminance of the adapting field, the whitepoint is found
    # at (100, 0, 0).
    L_A = np.inf
    cam16 = colorio.cs.CAM16UCS(0.69, 20, L_A)
    out = cam16.from_xyz100(colorio.illuminants.whitepoints_cie1931["D65"])
    assert np.all(out == [100, 0, 0])


if __name__ == "__main__":
    test_conversion(np.random.rand(3))
