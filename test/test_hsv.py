import numpy
import pytest

import colorio

numpy.random.seed(2)


@pytest.mark.parametrize(
    "vals,ref",
    [
        ([10, 20, 30], [210, 2 / 3, 2 / 17]),
        ([110, 120, 130], [210, 2 / 13, 26 / 51]),
        ([210, 220, 230], [210, 2 / 23, 46 / 51]),
    ],
)
def test_reference_srgb(vals, ref):
    hsv = colorio.Hsv()
    assert numpy.all(abs(hsv.from_srgb256(vals) - ref) < 1.0e-14 * numpy.array(ref))
    return


@pytest.mark.parametrize(
    "vals", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion(vals):
    hsv = colorio.Hsv()
    out = hsv.to_srgb1(hsv.from_srgb1(vals))
    assert numpy.all(abs(vals - out) < 1.0e-14)
    return
