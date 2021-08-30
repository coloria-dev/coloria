import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz,ref",
    [
        ([3.0, 4.0, 5.0], [3.0906319857540576, 3.9946283071224067, 5.3327338111342675]),
        ([67.3, 20.3, 71.9], [68.5805744415844, 20.85573438475875, 76.96296751820829]),
        (
            [90.0, 92.0, 96.0],
            [92.03403265014609, 92.01151985924933, 102.39979099626261],
        ),
    ],
)
def test_bfd(xyz, ref):
    cat, _ = colorio.cat.bradford(
        colorio.illuminants.whitepoints_cie1931["D65"],
        colorio.illuminants.whitepoints_cie1964["C"],
    )
    out = cat @ xyz

    print(list(out))
    assert np.all(np.abs(out - ref) < 1.0e-13 * out)
