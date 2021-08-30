import numpy as np

import colorio


def test_cmccat97():
    F = 1.0
    L_A = 20.0
    xyz = [67.3, 20.3, 71.9]
    whitepoint_test = colorio.illuminants.whitepoints_cie1931["D65"]
    whitepoint_reference = colorio.illuminants.whitepoints_cie1964["C"]
    ref = [68.57996639633065, 20.830994233615545, 77.36920346230376]

    cat = colorio.cat.CMCCAT97(whitepoint_test, whitepoint_reference, F, L_A)
    out = cat.apply(xyz)

    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)

    # check backwards transformation
    xyz2 = cat.apply_inv(out)
    assert np.all(np.abs(xyz - xyz2) < 1.0e-13 * out)
