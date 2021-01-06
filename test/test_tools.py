import numpy
import pytest

import colorio


def test_flat_gamut(xy_to_2d=lambda xy: xy):
    colorio.show_flat_gamut()


@pytest.mark.parametrize(
    "a", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion_variants(a):
    b = a + 1.0e-3 * numpy.random.rand(*a.shape)
    diff = colorio.delta(a, b)
    assert diff.shape == a.shape[1:]


def test_luo_rigg():
    colorio.show_luo_rigg(ellipse_scaling=1.5)


def test_xy_gamut_mesh():
    points, cells = colorio.xy_gamut_mesh(0.05)

    # import meshio
    # meshio.write_points_cells("test.vtu", points, {"triangle": cells})
    # exit(1)


if __name__ == "__main__":
    test_luo_rigg()
    # test_xy_gamut_mesh()
    # colorspace_ = colorio.SrgbLinear()
    # colorspace_ = colorio.Hdr()
    # colorspace_ = colorio.XYZ()
    # colorspace_ = colorio.XYY()
    # colorspace_ = colorio.IPT()
    # colorspace_ = colorio.JzAzBz()
    # colorspace_ = colorio.CIELUV()
    # colorspace_ = colorio.CIELAB()
    # colorspace_ = colorio.CAM02('UCS', 0.69, 20, 64/numpy.pi/5)
    # colorspace_ = colorio.CAM16UCS(0.69, 20, 64 / numpy.pi / 5)
    # test_hdr_gamut(colorspace_, n=10)
    # test_visible_gamut(colorspace_, cut_000=False)
    # test_srgb_gamut(colorspace_, cut_000=False)
