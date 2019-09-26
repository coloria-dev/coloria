import numpy
import pytest

import colorio


numpy.random.seed(0)


@pytest.mark.parametrize(
    "xyz", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion(xyz):
    rlab = colorio.RLAB()
    out = rlab.to_xyz100(rlab.from_xyz100(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([10, 20, 30], [49.56203373, -53.01687639, -13.37480413]),
        ([80, 90, 10], [95.27583974, -14.70510349, 101.4152633]),
        ([0.5, 0.6, 0.4], [10.78695859, -2.86443652, 3.39817955]),
    ],
)
def test_reference_xyz(xyz100, ref):
    rlab = colorio.RLAB()
    xyz100 = numpy.array(xyz100)
    assert numpy.all(
        abs(rlab.from_xyz100(xyz100) - ref) < 1.0e-4 * abs(numpy.array(ref))
    )
    return
