import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

import colorio


def test_flat_gamut():
    colorio.plot_xy_gamut().show()


def test_xy_gamut_mesh():
    # points, cells =
    colorio.xy_gamut_mesh(0.05)
    # import meshio
    # meshio.write_points_cells("test.vtu", points, {"triangle": cells})


@pytest.mark.parametrize(
    "cs,lightness",
    [
        [colorio.cs.XYY(1), 0.4],
        [colorio.cs.CIELAB(), 50],
        [colorio.cs.CAM16UCS(0.69, 20, 4.074), 50],
    ],
)
def test_visible_slice(cs, lightness):
    colorio.plot_visible_slice(cs, lightness)
    plt.show()
    # colorio.save_visible_slice("visible-slice.png", cs, lightness)


@pytest.mark.skip("Fails on gh-actions")
def test_rgb_slice():
    cs = colorio.cs.CIELAB()
    colorio.show_rgb_slice(cs, lightness=50, n=100)


@pytest.mark.parametrize(
    "colorspace",
    [
        colorio.cs.CIELAB(),
        # test against xyy to trigger ~self.is_origin_well_defined
        colorio.cs.XYY(1),
        colorio.cs.CAM02("UCS", 0.69, 20, 20),
    ],
)
@pytest.mark.parametrize("variant", ["srgb", "hdr"])
def test_srgb_gamut(colorspace, variant, n=10):
    with tempfile.TemporaryDirectory() as tmpdir:
        colorio.save_rgb_gamut(Path(tmpdir) / "srgb.vtu", colorspace, variant, n=n)


@pytest.mark.parametrize(
    "colorspace",
    [
        colorio.cs.CIELAB(),
        colorio.cs.XYY(1),
        colorio.cs.CAM02("UCS", 0.69, 20, 20),
    ],
)
def test_visible_gamut(colorspace):
    observer = colorio.observers.cie_1931_2()
    colorio.plot_visible_gamut(colorspace, observer, max_Y1=1)
    # show() segfaults when used in tests
    # p.show(screenshot=Path(tmpdir) / "cone.png")


def test_surface_gamut():
    colorspace = colorio.cs.XYY(1)
    illuminant = colorio.illuminants.d65()
    observer = colorio.observers.cie_1931_2()
    colorio.plot_surface_gamut(colorspace, observer, illuminant)


def test_srgb_gradient():
    # cs = colorio.cs.CIELAB()
    cs = colorio.cs.OKLAB()
    # cs = colorio.cs.DIN99()
    # cs = colorio.cs.CIELUV()
    colorio.plot_srgb255_gradient(cs, [0, 0, 255], [255, 255, 0]).show()


def test_primary_srgb_gradient():
    cs = colorio.cs.CIELAB()
    colorio.plot_primary_srgb_gradients(cs).show()
