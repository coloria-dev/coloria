import numpy as np

import colorio


def test_von_kries():
    xyz = [67.3, 20.3, 71.9]
    whitepoint_test = colorio.illuminants.whitepoints_cie1931["D65"]
    whitepoint_reference = colorio.illuminants.whitepoints_cie1964["C"]
    ref = [68.367018961615, 20.3, 76.69540240441576]

    cat = colorio.cat.VonKries(whitepoint_test, whitepoint_reference)
    out = cat.apply(xyz)

    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)

    # check backwards tranformation
    xyz2 = cat.apply_inv(out)
    assert np.all(np.abs(xyz - xyz2) < 1.0e-13 * out)
