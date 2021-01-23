import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        (
            [10.0, 20.0, 30.0],
            [0.5640791324061505, -0.17363647246705063, -0.03757807340002353],
        ),
        (
            [80.0, 90.0, 10.0],
            [0.9619070925241717, -0.04412900436798305, 0.2075039070964543],
        ),
        (
            [0.5, 0.6, 0.4],
            [0.1803246353234263, -0.010477821543242925, 0.013831761402855879],
        ),
        (
            colorio.illuminants.whitepoints_cie1931["D65"],
            [1.0, 0.0, 0.0],
        ),
    ],
)
def test_reference_xyz(xyz100, ref):
    cs = colorio.cs.OKLAB()
    xyz100 = np.asarray(xyz100)
    val = cs.from_xyz100(xyz100)
    print(ref)
    print(val.tolist())
    assert np.all(np.abs(val - ref) < 1.0e-4 * np.abs(ref) + 1.0e-4)


if __name__ == "__main__":
    cs = colorio.cs.OKLAB()
    # cs.save_srgb_gamut("oklab.vtk", n=50)
