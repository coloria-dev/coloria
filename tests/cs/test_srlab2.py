import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([10, 20, 30], [51.232467008097125, -47.92421786981609, -14.255014329381225]),
        ([80, 90, 10], [96.26425102758816, -30.92082858867076, 103.76703583290106]),
        (
            colorio.illuminants.whitepoints_cie1931["D65"],
            [99.99977248346777, -0.004069281557519844, -0.00039226988315022027],
        ),
        # make sure we have a test point with values below the linearity point
        ([0.5, 0.6, 0.4], [5.423523417804045, -2.6161648383355214, 3.6349016311770272]),
    ],
)
def test_reference_xyz(xyz100, ref):
    cs = colorio.cs.SRLAB2()
    vals = cs.from_xyz100(xyz100)
    print(list(vals))
    assert np.all(np.abs(vals - ref) < 1.0e-13 * np.abs(ref) + 1.0e-12)


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        (
            [2.61219, 1.52732, 10.96471],
            [13.596107145746988, -0.41689944789488464, -50.92072274345597],
        ),
        (
            [2.39318, 2.01643, 1.64315],
            [15.57318944741349, 9.84652234536864, 0.2249224382449384],
        ),
        (
            colorio.illuminants.whitepoints_cie1931["D50"],
            [99.99977248346777, -0.004069281557519844, -0.00039226988315022027],
        ),
        # make sure we have a test point with values below the linearity point
        ([0.5, 0.6, 0.4], [5.401958805705379, -3.1016033258418894, 1.75872187951029]),
    ],
)
def test_reference_xyz_d50(xyz100, ref):
    cs = colorio.cs.SRLAB2(whitepoint=colorio.illuminants.whitepoints_cie1931["D50"])
    vals = cs.from_xyz100(xyz100)
    print(list(vals))
    assert np.all(np.abs(vals - ref) < 1.0e-13 * np.abs(ref) + 1.0e-12)
