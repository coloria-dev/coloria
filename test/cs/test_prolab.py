import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([10, 20, 30], [64.92024825, -24.70816673, -7.55216015]),
        ([80, 90, 10], [98.1102344, 1.74587423, 37.37866513]),
        ([0.5, 0.6, 0.4], [3.77801345, -0.22163694, 0.59666459]),
        (
            colorio.illuminants.whitepoints_cie1931["D65"],
            [100.00063901, 4.72401578, 3.10710416],
        ),
        (
            colorio.illuminants.whitepoints_cie1931["D65"] / 100,
            [6.77847542, 0.3202142, 0.21061295],
        ),
    ],
)
def test_reference_xyz(xyz100, ref):
    cs = colorio.cs.PROLAB()
    xyz100 = np.asarray(xyz100)
    print(cs.from_xyz100(xyz100))
    assert np.all(np.abs(cs.from_xyz100(xyz100) - ref) < 1.0e-4 * np.abs(ref) + 1.0e-4)
