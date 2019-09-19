import numpy
import pytest

import colorio


@pytest.mark.parametrize(
    "xyz", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion(xyz):
    cielab = colorio.CIELAB()
    out = cielab.to_xyz100(cielab.from_xyz100(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([10, 20, 30], [51.8372, -56.3591, -13.1812]),
        ([80, 90, 10], [95.9968, -10.6593, 102.8625]),
        # make sure we have a test point with values below the linearity point
        ([0.5, 0.6, 0.4], [5.4198, -2.8790, 3.6230]),
    ],
)
def test_reference_xyz(xyz100, ref):
    cielab = colorio.CIELAB()
    xyz100 = numpy.array(xyz100)
    assert numpy.all(
        abs(cielab.from_xyz100(xyz100) - ref) < 1.0e-4 * abs(numpy.array(ref))
    )
    return


@pytest.mark.parametrize(
    "xyz100,ref",
    [
        ([2.61219, 1.52732, 10.96471], [12.7806, 26.1147, -52.4348]),
        ([2.39318, 2.01643, 1.64315], [15.5732, 9.7574, 0.2281]),
        # make sure we have a test point with values below the linearity point
        ([0.5, 0.6, 0.4], [5.4198, -3.1711, 1.7953]),
    ],
)
def test_reference_xyz_d50(xyz100, ref):
    cielab = colorio.CIELAB(whitepoint=colorio.illuminants.whitepoints_cie1931["D50"])
    xyz100 = numpy.array(xyz100)
    assert numpy.all(
        abs(cielab.from_xyz100(xyz100) - ref) < 1.0e-4 * abs(numpy.array(ref))
    )
    return


def test_macadams():
    # cielab = colorio.CIELAB()
    # cielab.show_macadams(0, 50)
    # cieluv = colorio.CIELUV()
    # cieluv.show_macadams(0, 50)
    jzazbz = colorio.JzAzBz()
    # jzazbz.show_macadams(0, 0.1)
    jzazbz.show_luo_rigg(0, 0.1)
    # xyy = colorio.XYY()
    # xyy.show_macadams(1.5, k0=2)
    #
    # L_A = 64 / numpy.pi / 5
    # cam02 = colorio.CAM02("UCS", 0.69, 20, L_A)
    # cam02.show_macadams(0, 50)
    #
    # L_A = 64 / numpy.pi / 5
    # cam16 = colorio.CAM16UCS(0.69, 20, L_A)
    # # cam16.show_macadams(0, 60)
    # cam16.show_luo_rigg(0, 60, ellipse_scaling=2.0)
    return


if __name__ == "__main__":
    test_macadams()
