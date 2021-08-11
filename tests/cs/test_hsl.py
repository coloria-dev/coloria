import numpy as np
import pytest

import colorio

rng = np.random.default_rng(2)


@pytest.mark.parametrize(
    "vals,ref",
    [
        ([10, 20, 30], [210, 0.5, 4 / 51]),
        ([110, 120, 130], [210, 25 / 300, 8 / 17]),
        ([210, 220, 230], [210, 2 / 7, 44 / 51]),
    ],
)
def test_reference_srgb(vals, ref):
    hsl = colorio.cs.HSL()
    assert np.all(abs(hsl.from_rgb256(vals) - ref) < 1.0e-14 * np.array(ref))


@pytest.mark.parametrize(
    "vals",
    [
        rng.random(3),
        rng.random((3, 7)),
        rng.random((3, 4, 5)),
        [1.0, 1.0, 1.0],
    ],
)
def test_conversion(vals):
    hsl = colorio.cs.HSL()
    out = hsl.to_rgb1(hsl.from_rgb1(vals))
    assert np.all(abs(vals - out) < 1.0e-14)


def test_nan():
    hsl = colorio.cs.HSL()

    vals = np.full(3, np.nan)
    out = hsl.from_rgb1(vals)
    print(out)
    assert np.all(np.isnan(out))

    vals = np.full(3, np.nan)
    out = hsl.to_rgb1(vals)
    assert np.all(np.isnan(out))
