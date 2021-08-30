import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz,ref",
    [
        ([3.0, 4.0, 5.0], [3.0822082672470015, 3.997846177980249, 5.285612492854658]),
        (
            [67.3, 20.3, 71.9],
            [68.20181638730817, 20.75223632965131, 75.97638883814696],
        ),
        (
            [90.0, 92.0, 96.0],
            [91.73039347580564, 91.99912497671738, 101.50017886075341],
        ),
    ],
)
def test_reference_value(xyz, ref):
    cat, _ = colorio.cat.cat02(
        colorio.illuminants.whitepoints_cie1931["D65"],
        colorio.illuminants.whitepoints_cie1964["C"],
        F=1.0,
        L_A=20.0,
    )
    out = cat @ xyz
    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)
