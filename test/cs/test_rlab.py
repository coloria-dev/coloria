import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([10, 20, 30], [49.56203373, -53.01687639, -13.37480413]),
        ([80, 90, 10], [95.27583974, -14.70510349, 101.4152633]),
        ([0.5, 0.6, 0.4], [10.78695859, -2.86443652, 3.39817955]),
    ],
)
def test_reference_xyz(xyz100, ref):
    rlab = colorio.cs.RLAB()
    xyz100 = np.array(xyz100)
    assert np.all(np.abs(rlab.from_xyz100(xyz100) - ref) < 1.0e-4 * np.abs(ref))
