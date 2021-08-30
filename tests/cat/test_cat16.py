import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz,ref",
    [
        ([3.0, 4.0, 5.0], [3.0792732980061093, 3.998788227557918, 5.284802632938603]),
        ([67.3, 20.3, 71.9], [68.20137123414864, 20.318197018662612, 75.8855939879793]),
        (
            [90.0, 92.0, 96.0],
            [91.74397001307572, 92.00516520876931, 101.50609654477812],
        ),
    ],
)
def test_reference_value(xyz, ref):
    cat, _ = colorio.cat.cat16(
        colorio.illuminants.whitepoints_cie1931["D65"],
        colorio.illuminants.whitepoints_cie1964["C"],
        F=1.0,
        L_A=20.0,
    )
    out = cat @ xyz
    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)
