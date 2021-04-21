import numpy as np
import pytest

import colorio

rng = np.random.default_rng(2)


@pytest.mark.parametrize(
    "vals,ref",
    [
        ([10, 20, 30], [210, 2 / 3, 2 / 17]),
        ([110, 120, 130], [210, 2 / 13, 26 / 51]),
        ([210, 220, 230], [210, 2 / 23, 46 / 51]),
    ],
)
def test_reference_srgb(vals, ref):
    hsv = colorio.cs.HSV()
    assert np.all(abs(hsv.from_srgb256(vals) - ref) < 1.0e-14 * np.array(ref))


@pytest.mark.parametrize(
    "vals", [rng.random(3), rng.random((3, 7)), rng.random((3, 4, 5))]
)
def test_conversion(vals):
    hsv = colorio.cs.HSV()
    out = hsv.to_srgb1(hsv.from_srgb1(vals))
    assert np.all(abs(vals - out) < 1.0e-14)
