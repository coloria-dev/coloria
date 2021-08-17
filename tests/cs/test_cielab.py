import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([10, 20, 30], [51.83721152653848, -56.35909588563908, -13.181163974878515]),
        ([80, 90, 10], [95.99676861425304, -10.659336431519876, 102.86254166243415]),
        (colorio.illuminants.whitepoints_cie1931["D65"], [100.0, 0.0, 0.0]),
        # make sure we have a test point with values below the linearity point
        ([0.5, 0.6, 0.4], [5.419777777777778, -2.87904161644785, 3.623046586533081]),
    ],
)
def test_reference_xyz(xyz100, ref):
    cs = colorio.cs.CIELAB()
    vals = cs.from_xyz100(xyz100)
    print(list(vals))
    assert np.all(np.abs(vals - ref) < 1.0e-12 * np.abs(ref) + 1.0e-12)


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        (
            [2.61219, 1.52732, 10.96471],
            [12.780698868103496, 26.11400692261134, -52.43465583230822],
        ),
        (
            [2.39318, 2.01643, 1.64315],
            [15.573231888842574, 9.757357302738868, 0.2280840026762867],
        ),
        (colorio.illuminants.whitepoints_cie1931["D50"], [100.0, 0.0, 0.0]),
        # make sure we have a test point with values below the linearity point
        ([0.5, 0.6, 0.4], [5.419777777777778, -3.1711206601843607, 1.7952998271595977]),
    ],
)
def test_reference_xyz_d50(xyz100, ref):
    cs = colorio.cs.CIELAB(whitepoint=colorio.illuminants.whitepoints_cie1931["D50"])
    vals = cs.from_xyz100(xyz100)
    print(list(vals))
    assert np.all(np.abs(vals - ref) < 1.0e-12 * np.abs(ref) + 1.0e-12)
