import numpy as np
import pytest
from cam16_legacy import CAM16Legacy

np.random.seed(0)


@pytest.mark.parametrize(
    "xyz",
    [
        100 * np.random.rand(3),
        100 * np.random.rand(3, 7),
        100 * np.random.rand(3, 4, 5),
    ],
)
def test_conversion(xyz):
    # test with srgb conditions
    L_A = 64 / np.pi / 5
    cam16 = CAM16Legacy(0.69, 20, L_A)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    out = cam16.to_xyz100(np.array([J, C, H]), "JCH")
    assert np.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = cam16.to_xyz100(np.array([Q, M, h]), "QMh")
    assert np.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = cam16.to_xyz100(np.array([J, s, h]), "Jsh")
    assert np.all(abs(xyz - out) < 1.0e-13 * abs(xyz))


if __name__ == "__main__":
    test_conversion(np.random.rand(3))
