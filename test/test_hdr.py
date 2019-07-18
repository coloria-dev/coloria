import numpy
import pytest

import colorio


@pytest.mark.parametrize(
    "vals", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion(vals):
    cs = colorio.Hdr()

    out = cs.to_xyz100(cs.from_xyz100(vals))
    assert numpy.all(abs(vals - out) < 1.0e-14)

    out = cs.to_hdr1(cs.from_hdr1(vals))
    assert numpy.all(abs(vals - out) < 1.0e-14)
    return
