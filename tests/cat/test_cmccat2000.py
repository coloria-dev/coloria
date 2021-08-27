import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz,ref",
    [
        (
            [3.0, 4.0, 5.0],
            [3.0873052916804333, 3.99630424803264, 5.319870094116135],
        ),
        (
            [67.3, 20.3, 71.9],
            [68.44010878974217, 20.649791873827574, 76.43399830620979],
        ),
        (
            [90.0, 92.0, 96.0],
            [91.95798638491023, 92.00870538929831, 102.16434737446967],
        ),
    ],
)
def test_reference_value(xyz, ref):
    cat = colorio.cat.CMCCAT2000(
        F=1.0,
        L_A1=20.0,
        L_A2=30.0,
        whitepoint_test=colorio.illuminants.whitepoints_cie1931["D65"],
        whitepoint_reference=colorio.illuminants.whitepoints_cie1964["C"],
    )
    out = cat.apply(xyz)
    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)
