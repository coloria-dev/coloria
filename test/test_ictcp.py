import numpy
import pytest

import colorio


@pytest.mark.parametrize(
    "rgb", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion(rgb):
    cs = colorio.ICtCp()
    out = cs.to_rec2100(cs.from_rec2100(rgb))
    assert numpy.all(abs(rgb - out) < 1.0e-8 * rgb)
    return
