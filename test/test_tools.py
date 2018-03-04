# -*- coding: utf-8 -*-
#
import numpy
import pytest

import colorio


@pytest.mark.parametrize('colorspace, cut_000', [
    # colorio.CIELAB(),
    (colorio.XYY(), True),
    (colorio.CAM02('UCS', 0.69, 20, 64/numpy.pi/5), False),
    ])
def test_srgb_gamut(colorspace, cut_000, n=10):
    colorio.show_srgb_gamut(colorspace, 'srgb.vtu', n=n, cut_000=cut_000)
    return


@pytest.mark.parametrize('colorspace', [
    colorio.CIELAB(),
    colorio.CAM02('UCS', 0.69, 20, 64/numpy.pi/5),
    ])
def test_hdr_gamut(colorspace, n=10):
    colorio.show_hdr_gamut(colorspace, 'hdr.vtu', n=n)
    return


@pytest.mark.parametrize('colorspace,cut_000', [
    # (colorio.CIELAB(), False),
    (colorio.XYY(), True),
    (colorio.CAM02('UCS', 0.69, 20, 64/numpy.pi/5), False),
    ])
def test_visible_gamut(colorspace, cut_000):
    illuminant = colorio.illuminants.d65()
    observer = colorio.observers.cie_1931_2()
    colorio.show_visible_gamut(
        colorspace, observer, illuminant, 'visible.vtu', cut_000=cut_000
        )
    return


def test_gamut_diagram():
    colorio.show_gamut_diagram()
    return


@pytest.mark.parametrize('a', [
    numpy.random.rand(3),
    numpy.random.rand(3, 7),
    numpy.random.rand(3, 4, 5),
    ])
def test_conversion_variants(a):
    b = a + 1.0e-3 * numpy.random.rand(*a.shape)
    diff = colorio.delta(a, b)
    assert diff.shape == a.shape[1:]
    return


@pytest.mark.parametrize('colorspace', [
    colorio.CIELAB(),
    colorio.XYY(),
    colorio.CAM02('UCS', 0.69, 20, 64/numpy.pi/5),
    ])
def test_ebner_fairchild(colorspace):
    colorio.show_ebner_fairchild(colorspace)
    return


@pytest.mark.parametrize('colorspace', [
    colorio.CIELAB(),
    colorio.XYY(),
    colorio.CAM02('UCS', 0.69, 20, 64/numpy.pi/5),
    ])
def test_hung_berns(colorspace):
    colorio.show_hung_berns(colorspace)
    return


@pytest.mark.parametrize('colorspace', [
    colorio.CIELAB(),
    ])
def test_munsell(colorspace):
    colorio.show_munsell(colorspace, V=5)
    return


def test_macadam():
    colorio.show_macadam(scaling=10)
    return


def test_luo_rigg():
    colorio.show_luo_rigg(scaling=1.5)
    return


if __name__ == '__main__':
    # test_luo_rigg()
    test_macadam()
    exit(1)
    # colorspace_ = colorio.SrgbLinear()
    # colorspace_ = colorio.Rec2020()
    # colorspace_ = colorio.XYZ()
    # colorspace_ = colorio.XYY()
    # colorspace_ = colorio.IPT()
    colorspace_ = colorio.JzAzBz()
    # colorspace_ = colorio.CIELUV()
    # colorspace_ = colorio.CIELAB()
    # colorspace_ = colorio.CAM02('UCS', 0.69, 20, 64/numpy.pi/5)
    # colorspace_ = colorio.CAM16UCS(0.69, 20, 64/numpy.pi/5)
    # test_hdr_gamut(colorspace_, n=10)
    # test_visible_gamut(colorspace_, cut_000=False)
    # test_srgb_gamut(colorspace_, cut_000=False)
    # test_ebner_fairchild(colorspace_)
    test_hung_berns(colorspace_)
    # test_munsell(colorspace_)
