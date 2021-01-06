import numpy
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion(xyz):
    cielab = colorio.CIELAB()
    out = cielab.to_xyz100(cielab.from_xyz100(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([10, 20, 30], [51.8372, -56.3591, -13.1812]),
        ([80, 90, 10], [95.9968, -10.6593, 102.8625]),
        (colorio.illuminants.whitepoints_cie1931["D65"], [100.0, 0.0, 0.0]),
        # make sure we have a test point with values below the linearity point
        ([0.5, 0.6, 0.4], [5.4198, -2.8790, 3.6230]),
    ],
)
def test_reference_xyz(xyz100, ref):
    cielab = colorio.CIELAB()
    xyz100 = numpy.array(xyz100)
    assert numpy.all(
        numpy.abs(cielab.from_xyz100(xyz100) - ref) < 1.0e-4 * numpy.abs(ref) + 1.0e-15
    )


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([2.61219, 1.52732, 10.96471], [12.7806, 26.1147, -52.4348]),
        ([2.39318, 2.01643, 1.64315], [15.5732, 9.7574, 0.2281]),
        (colorio.illuminants.whitepoints_cie1931["D50"], [100.0, 0.0, 0.0]),
        # make sure we have a test point with values below the linearity point
        ([0.5, 0.6, 0.4], [5.4198, -3.1711, 1.7953]),
    ],
)
def test_reference_xyz_d50(xyz100, ref):
    cielab = colorio.CIELAB(whitepoint=colorio.illuminants.whitepoints_cie1931["D50"])
    xyz100 = numpy.array(xyz100)
    assert numpy.all(
        numpy.abs(cielab.from_xyz100(xyz100) - ref) < 1.0e-4 * numpy.abs(ref) + 1.0e-15
    )
