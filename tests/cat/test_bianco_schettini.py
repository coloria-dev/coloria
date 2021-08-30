import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz,ref",
    [
        ([3.0, 4.0, 5.0], [3.091013955314862, 3.993389626548669, 5.333428529739702]),
        (
            [67.3, 20.3, 71.9],
            [68.45572007411377, 20.995067476918063, 76.63085190939026],
        ),
        (
            [90.0, 92.0, 96.0],
            [92.03525718564795, 92.01390843134203, 102.40453352447496],
        ),
    ],
)
def test_bianco_schettini(xyz, ref):
    cat, _ = colorio.cat.bianco_schettini(
        colorio.illuminants.whitepoints_cie1931["D65"],
        colorio.illuminants.whitepoints_cie1964["C"],
    )
    out = cat @ xyz

    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)


@pytest.mark.parametrize(
    "xyz,ref",
    [
        ([3.0, 4.0, 5.0], [3.0847412300813093, 3.9954689664563747, 5.332350480399397]),
        (
            [67.3, 20.3, 71.9],
            [68.44211847845229, 20.366887647522212, 76.40174168110818],
        ),
        ([90.0, 92.0, 96.0], [92.06457641666026, 92.0193412804294, 102.41500324735371]),
    ],
)
def test_bianco_schettini_pos(xyz, ref):
    cat, _ = colorio.cat.bianco_schettini_pos(
        colorio.illuminants.whitepoints_cie1931["D65"],
        colorio.illuminants.whitepoints_cie1964["C"],
    )
    out = cat @ xyz

    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)
