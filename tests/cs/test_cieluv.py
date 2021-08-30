import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([10, 20, 30], [51.83721152653848, -65.9326680348282, -12.356536578386144]),
        ([80, 90, 10], [95.99676861425304, 26.62924882039604, 107.89622364611081]),
        (colorio.illuminants.whitepoints_cie1931["D65"], [100.0, 0.0, 0.0]),
        # make sure we have a test point with values below the linearity point
        ([0.5, 0.6, 0.4], [5.419777777777778, -0.7696690409323268, 2.560171459447575]),
    ],
)
def test_reference_xyz(xyz100, ref):
    cs = colorio.cs.CIELUV()
    vals = cs.from_xyz100(xyz100)
    print(xyz100)
    print(list(vals))
    assert np.all(np.abs(vals - ref) < 1.0e-12 * np.abs(ref) + 1.0e-12)


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        (
            [2.61219, 1.52732, 10.96471],
            [12.780698868103496, -5.033038257263479, -41.99645404861298],
        ),
        (
            [2.39318, 2.01643, 1.64315],
            [15.573231888842574, 9.240537370066974, -1.0163158395602911],
        ),
        (colorio.illuminants.whitepoints_cie1931["D50"], [100.0, 0.0, 0.0]),
        # make sure we have a test point with values below the linearity point
        ([0.5, 0.6, 0.4], [5.419777777777778, -1.5672596028594659, 1.169553707268332]),
    ],
)
def test_reference_xyz_d50(xyz100, ref):
    cs = colorio.cs.CIELUV(whitepoint=colorio.illuminants.whitepoints_cie1931["D50"])
    vals = cs.from_xyz100(xyz100)
    print(xyz100)
    print(list(vals))
    assert np.all(np.abs(vals - ref) < 1.0e-12 * np.abs(ref) + 1.0e-12)
