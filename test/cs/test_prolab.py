import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz1,ref",
    [
        # Reference values from
        # https://github.com/konovalenko-iitp/proLab/issues/1
        ([0, 0, 0], [0, 0, 0]),
        ([0.4125, 0.2127, 0.0193], [63.8813, 64.6958, 26.6611]),
        ([0.3576, 0.7152, 0.1192], [93.2702, -46.2794, 30.9420]),
        ([0.1804, 0.0722, 0.9503], [67.8217, 19.0252, -65.5766]),
        ([0.5380, 0.7873, 1.0695], [96.5227, -23.0585, -9.5024]),
        ([0.5929, 0.2848, 0.9696], [82.3994, 48.2089, -32.2321]),
        ([0.7700, 0.9278, 0.1385], [98.6665, -10.1353, 34.9347]),
        ([0.9505, 1.0000, 1.0888], [100.0000, -0.0000, -0.0000]),
    ],
)
def test_reference_xyz(xyz1, ref):
    xyz100 = np.array(xyz1) * 100
    cs = colorio.cs.PROLAB()
    xyz100 = np.asarray(xyz100)
    assert np.all(np.abs(cs.from_xyz100(xyz100) - ref) < 1.0e-4 * np.abs(ref) + 1.0e-4)
