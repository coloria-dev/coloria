import numpy
import pytest

from cam16_legacy import CAM16Legacy

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
    cam16 = CAM16Legacy(0.69, 20, L_A)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    out = cam16.to_xyz100(numpy.array([J, C, H]), "JCH")
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = cam16.to_xyz100(numpy.array([Q, M, h]), "QMh")
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = cam16.to_xyz100(numpy.array([J, s, h]), "Jsh")
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))
    return


if __name__ == "__main__":
    test_conversion(numpy.random.rand(3))
