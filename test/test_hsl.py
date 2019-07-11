import numpy
import pytest

import colorio

# @pytest.mark.parametrize(
#     "vals", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
# )
# def test_conversion(vals):
#     hsl = colorio.Hsl()
#     out = hsl.to_srgb1(hsl.from_srgb1(vals))
#     assert numpy.all(abs(vals - out) < 1.0e-14)
#     return


@pytest.mark.parametrize(
    "vals,ref",
    [
        ([10, 20, 30], [210.0, 50.0, 400 / 51]),
        ([110, 120, 130], [210.0, 25 / 3, 47.05882352941176]),
        ([210, 220, 230], [210.0, 200 / 7, 4400 / 51]),
    ],
)
def test_reference_srgb(vals, ref):
    hsl = colorio.Hsl()
    assert numpy.all(abs(hsl.from_srgb256(vals) - ref) < 1.0e-14 * numpy.array(ref))
    return
