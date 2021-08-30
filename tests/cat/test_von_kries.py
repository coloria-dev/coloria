import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz,ref",
    [
        ([3.0, 4.0, 5.0], [3.1016365089661244, 4.0, 5.333477218665907]),
        (
            [67.3, 20.3, 71.9],
            [68.367018961615, 20.3, 76.69540240441576],
        ),
        ([90.0, 92.0, 96.0], [91.9882823645901, 92.0, 102.40276259838542]),
    ],
)
def test_von_kries(xyz, ref):
    cat, _ = colorio.cat.von_kries(
        colorio.illuminants.whitepoints_cie1931["D65"],
        colorio.illuminants.whitepoints_cie1964["C"],
    )
    out = cat @ xyz

    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)
