import numpy
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([10, 20, 30], [2.6182234, -0.80594911, -0.1744219]),
        ([80, 90, 10], [4.46477722, -0.20482869, 0.96314782]),
        ([0.5, 0.6, 0.4], [0.83699281, -0.04863374, 0.06420135]),
        (
            colorio.illuminants.whitepoints_cie1931["D65"],
            [4.64158795e00, -4.68417772e-05, -3.99681552e-04],
        ),
        (
            colorio.illuminants.whitepoints_cie1931["D65"] / 100,
            [1.0, 0.0, 0.0],
        ),
    ],
)
def test_reference_xyz(xyz100, ref):
    cs = colorio.cs.OKLAB()
    xyz100 = numpy.asarray(xyz100)
    assert numpy.all(
        numpy.abs(cs.from_xyz100(xyz100) - ref) < 1.0e-4 * numpy.abs(ref) + 1.0e-4
    )


if __name__ == "__main__":
    cs = colorio.cs.OKLAB()
    # cs.save_srgb_gamut("oklab.vtk", n=50)
    cs.show_ebner_fairchild()
