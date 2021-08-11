import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "vals",
    [
        100 * np.random.rand(3),
        100 * np.random.rand(3, 7),
        100 * np.random.rand(3, 4, 5),
    ],
)
def test_conversion(vals):
    cs = colorio.cs.HdrLinear()

    out = cs.to_rgb1(cs.from_rgb1(vals))
    assert np.all(abs(vals - out) < 1.0e-14)
