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
    cat, _ = colorio.cat.cmccat2000(
        colorio.illuminants.whitepoints_cie1931["D65"],
        colorio.illuminants.whitepoints_cie1964["C"],
        F=1.0,
        L_A1=20.0,
        L_A2=30.0,
    )
    out = cat @ xyz
    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)


def test_reference():
    # numbers from the article
    #
    # CMC 2000 Chromatic Adaptation Transform: CMCCAT2000
    # Changjun Li, M. Ronnier Luo,* Bryan Rigg, Robert W. G. Hunt
    #
    cat, _ = colorio.cat.cmccat2000(
        [111.15, 100.00, 35.20],
        [94.81, 100.00, 107.30],
        F=1.0,
        L_A1=200.0,
        L_A2=200.0,
    )
    out = cat @ [22.48, 22.74, 8.54]

    assert np.all(out.round(2) == [19.53, 23.07, 24.97])
