import numpy
import pytest

import colorio

numpy.random.seed(0)


@pytest.mark.parametrize(
    "xyz",
    [
        100 * numpy.random.rand(3),
        100 * numpy.random.rand(3, 7),
        100 * numpy.random.rand(3, 4, 5),
    ],
)
def test_conversion(xyz):
    # test with srgb conditions
    L_A = 64 / numpy.pi / 5
    cam16 = colorio.CAM16(0.69, 20, L_A)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    out = cam16.to_xyz100(numpy.array([J, C, H]), "JCH")
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = cam16.to_xyz100(numpy.array([Q, M, h]), "QMh")
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = cam16.to_xyz100(numpy.array([J, s, h]), "Jsh")
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))


@pytest.mark.parametrize("xyz", [numpy.zeros(3), numpy.zeros((3, 4, 5))])
def test_zero(xyz):
    cam16 = colorio.CAM16(0.69, 20, 64 / numpy.pi / 5)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    assert numpy.all(J == 0.0)
    assert numpy.all(C == 0.0)
    assert numpy.all(h == 0.0)
    assert numpy.all(M == 0.0)
    assert numpy.all(s == 0.0)
    assert numpy.all(Q == 0.0)

    out = cam16.to_xyz100(numpy.array([J, C, H]), "JCH")
    assert numpy.all(abs(out) < 1.0e-13)

    out = cam16.to_xyz100(numpy.array([Q, M, h]), "QMh")
    assert numpy.all(abs(out) < 1.0e-13)

    out = cam16.to_xyz100(numpy.array([J, s, h]), "Jsh")
    assert numpy.all(abs(out) < 1.0e-13)


@pytest.mark.parametrize(
    "xyz", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion_variants(xyz):
    # test with srgb conditions
    L_A = 64 / numpy.pi / 5
    cam16 = colorio.CAM16UCS(0.69, 20, L_A)
    out = cam16.to_xyz100(cam16.from_xyz100(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)


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
    L_A = 64 / numpy.pi / 5
    cam16 = colorio.CAM16UCS(0.69, 20, L_A)
    out = cam16.from_xyz100(xyz)
    ref = numpy.array(ref)
    assert numpy.all(abs(ref - out) < 1.0e-14 * ref)
    return


def test_whitepoint():
    # With infinite luminance of the adapting field, the whitepoint is found
    # at (100, 0, 0).
    L_A = numpy.inf
    cam16 = colorio.CAM16UCS(0.69, 20, L_A)
    out = cam16.from_xyz100(colorio.illuminants.whitepoints_cie1931["D65"])
    assert numpy.all(out == [100, 0, 0])
    return


if __name__ == "__main__":
    test_conversion(numpy.random.rand(3))
