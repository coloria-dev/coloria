import tempfile
from pathlib import Path

import numpy
import pytest

import colorio


@pytest.mark.parametrize(
    "cs,k0,level",
    [
        [colorio.cs.XYY(1), 2, 0.4],
        [colorio.cs.CIELAB(), 0, 50],
        [colorio.cs.CAM16UCS(0.69, 20, 4.074), 0, 50],
    ],
)
def test_visible_slice(cs, k0, level):
    cs.show_visible_slice(k0, level)
    # cs.save_visible_slice("visible-slice.png", k0, level)


@pytest.mark.parametrize("variant", ["srgb", "hdr"])
@pytest.mark.parametrize(
    "colorspace",
    [
        colorio.cs.CIELAB(),
        # colorio.cs.XYY(1),
        colorio.cs.CAM02("UCS", 0.69, 20, 64 / numpy.pi / 5),
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
        colorio.cs.CAM02("UCS", 0.69, 20, 64 / numpy.pi / 5),
    ],
)
def test_cone_gamut(colorspace, n=10):
    observer = colorio.observers.cie_1931_2()
    with tempfile.TemporaryDirectory() as tmpdir:
        colorspace.save_cone_gamut(Path(tmpdir) / "cone.vtu", observer, max_Y=1)


@pytest.mark.parametrize(
    "colorspace",
    [
        # colorio.cs.CIELAB(),
        # colorio.cs.XYY(1),
        colorio.cs.CAM02("UCS", 0.69, 20, 64 / numpy.pi / 5),
    ],
)
def test_visible_gamut(colorspace):
    illuminant = colorio.illuminants.d65()
    observer = colorio.observers.cie_1931_2()
    with tempfile.TemporaryDirectory() as tmpdir:
        colorspace.save_visible_gamut(
            observer, illuminant, Path(tmpdir) / "visible.vtu"
        )
