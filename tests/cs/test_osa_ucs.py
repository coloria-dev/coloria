import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz",
    [
        [0.0, 0.0, 0.0],
        # difficult case that fails if the initial values aren't chosen carefully
        [12.0, 67.0, 20.0],
    ],
)
def test_conversion(xyz):
    osa = colorio.cs.OsaUcs()
    out = osa.to_xyz100(osa.from_xyz100(xyz))
    assert np.all(abs(xyz - out) < 1.0e-9)
