import numpy as np

import colorio


def test_refereece_value():
    F = 1.0
    L_A1 = 20.0
    L_A2 = 30.0
    xyz = [67.3, 20.3, 71.9]
    whitepoint_test = colorio.illuminants.whitepoints_cie1931["D65"]
    whitepoint_reference = colorio.illuminants.whitepoints_cie1964["C"]
    ref = [68.57996639633065, 20.830994233615545, 77.36920346230376]

    cat = colorio.cat.CMCCAT2000(F, L_A1, L_A2, whitepoint_test, whitepoint_reference)
    out = cat.apply(xyz)
    print(list(out))

    assert np.all(np.abs(out - ref) < 1.0e-13 * out)
