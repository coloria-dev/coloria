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


def test_macadam():
    def xy_to_2d(xy):
        x, y = xy
        return numpy.array([4 * x, 9 * y]) / (-2 * x + 12 * y + 3)

    # def xy_to_2d(xy):
    #     return xy

    colorio.show_macadam(
        ellipse_scaling=10,
        xy_to_2d=xy_to_2d,
        # plot_standard_deviations=True,
        # axes_labels=['u\'', 'v\'']
    )


def test_luo_rigg():
    colorio.show_luo_rigg(ellipse_scaling=1.5)


def test_xy_gamut_mesh():
    points, cells = colorio.xy_gamut_mesh(0.05)

    # import meshio
    # meshio.write_points_cells("test.vtu", points, {"triangle": cells})
    # exit(1)


if __name__ == "__main__":
    # test_luo_rigg()
    # test_xy_gamut_mesh()
    # test_macadam()
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
    test_luo_rigg()
