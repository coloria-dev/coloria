import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz,ref",
    [
        ([3.0, 4.0, 5.0], [3.083180036109638, 3.9902781600694097, 5.332277989046864]),
        (
            [67.3, 20.3, 71.9],
            [68.75286586870737, 21.024453944792988, 76.82507941523572],
        ),
        ([90.0, 92.0, 96.0], [92.0643506038423, 92.0275863474639, 102.40520006943595]),
    ],
)
def test_bfd(xyz, ref):
    cat, _ = colorio.cat.sharp(
        colorio.illuminants.whitepoints_cie1931["D65"],
        colorio.illuminants.whitepoints_cie1964["C"],
    )
    out = cat @ xyz

    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)
