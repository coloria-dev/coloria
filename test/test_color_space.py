import pytest

import colorio


@pytest.mark.parametrize(
    "cs,k0,level",
    [
        [colorio.XYY(), 2, 0.4],
        [colorio.CIELAB(), 0, 50],
        [colorio.CAM16UCS(0.69, 20, 4.074), 0, 50],
    ],
)
def test_visible_slice(cs, k0, level):
    cs.show_visible_slice(k0, level)
    # cs.save_visible_slice("visible-slice.png", k0, level)


@pytest.mark.parametrize(
    "cs,k0,level",
    [[colorio.XYY(), 2, 0.4], [colorio.CIELUV(), 0, 50], [colorio.JzAzBz(), 0, 0.5]],
)
def test_luo_rigg(cs, k0, level):
    cs.show_luo_rigg(k0, level)
    cs.save_luo_rigg("luo-rigg.png", k0, level)
