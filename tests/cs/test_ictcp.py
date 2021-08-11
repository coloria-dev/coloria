import numpy as np
import pytest

import colorio

rng = np.random.default_rng(0)


@pytest.mark.parametrize(
    "rgb", [rng.random(3), rng.random((3, 7)), rng.random((3, 4, 5))]
)
def test_conversion(rgb):
    cs = colorio.cs.ICtCp()
    out = cs.to_rec2100(cs.from_rec2100(rgb))
    assert np.all(abs(rgb - out) < 1.0e-8 * rgb)


@pytest.mark.parametrize(
    "rgb",
    [
        100 * rng.random(3),
        100 * rng.random((3, 7)),
        100 * rng.random((3, 4, 5)),
    ],
)
def test_conversion_xyz100(rgb):
    xyz100 = colorio.cs.HdrLinear().to_xyz100(rgb)
    cs = colorio.cs.ICtCp()
    out = cs.to_xyz100(cs.from_xyz100(xyz100))
    assert np.all(abs(xyz100 - out) < 1.0e-8 * xyz100)
