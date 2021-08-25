import numpy as np
import pytest
from cam16_legacy import CAM16Legacy

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
    cam16 = CAM16Legacy(0.69, 20, 10)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    out = cam16.to_xyz100(np.array([J, C, H]), "JCH")
    assert np.all(abs(xyz - out) < tol * abs(xyz))

    out = cam16.to_xyz100(np.array([Q, M, h]), "QMh")
    assert np.all(abs(xyz - out) < tol * abs(xyz))

    out = cam16.to_xyz100(np.array([J, s, h]), "Jsh")
    assert np.all(abs(xyz - out) < tol * abs(xyz))


if __name__ == "__main__":
    test_conversion(rng.random(3))
