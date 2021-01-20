import tempfile
from pathlib import Path

import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "cs,lightness",
    [
        [colorio.cs.XYY(1), 0.4],
        [colorio.cs.CIELAB(), 50],
        [colorio.cs.CAM16UCS(0.69, 20, 4.074), 50],
    ],
)
def test_visible_slice(cs, lightness):
    cs.show_visible_slice(lightness)
    # cs.save_visible_slice("visible-slice.png", lightness)


def test_rgb_slice():
    cs = colorio.cs.CIELAB()
    lightness = 50
    cs.show_rgb_slice(lightness, n=100)


@pytest.mark.parametrize("variant", ["srgb", "hdr"])
@pytest.mark.parametrize(
    "colorspace",
    [
        colorio.cs.CIELAB(),
        # test against xyy to trigger ~self.is_origin_well_defined
        colorio.cs.XYY(1),
        colorio.cs.CAM02("UCS", 0.69, 20, 64 / np.pi / 5),
    ],
)
def test_srgb_gamut(variant, colorspace, n=10):
    with tempfile.TemporaryDirectory() as tmpdir:
        colorspace.save_rgb_gamut(Path(tmpdir) / "srgb.vtu", variant, n=n)


@pytest.mark.parametrize(
    "colorspace",
    [
        colorio.cs.CIELAB(),
        colorio.cs.XYY(1),
        colorio.cs.CAM02("UCS", 0.69, 20, 64 / np.pi / 5),
    ],
)
def test_cone_gamut(colorspace):
    observer = colorio.observers.cie_1931_2()
    with tempfile.TemporaryDirectory() as tmpdir:
        colorspace.save_cone_gamut(Path(tmpdir) / "cone.vtu", observer, max_Y=1)


def test_visible_gamut():
    colorspace = colorio.cs.XYY(1)
    illuminant = colorio.illuminants.d65()
    observer = colorio.observers.cie_1931_2()
    with tempfile.TemporaryDirectory() as tmpdir:
        colorspace.save_visible_gamut(
            observer, illuminant, Path(tmpdir) / "visible.vtu"
        )


def test_srgb_gradient():
    # cs = colorio.cs.CIELAB()
    cs = colorio.cs.OKLAB()
    # cs = colorio.cs.DIN99()
    # cs = colorio.cs.CIELUV()
    cs.show_srgb_gradient([0, 0, 255], [255, 255, 0])


def test_primary_srgb_gradient():
    cs = colorio.cs.CIELAB()
    # cs = colorio.cs.OKLAB()
    cs.show_primary_srgb_gradients()


if __name__ == "__main__":
    # test_srgb_gradient()
    test_primary_srgb_gradient()
