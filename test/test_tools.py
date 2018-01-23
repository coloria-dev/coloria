# -*- coding: utf-8 -*-
#
import numpy
import pytest

import colorio


@pytest.mark.parametrize('colorspace', [
    colorio.CIELAB(),
    colorio.CAM02('UCS', 0.69, 20, 64/numpy.pi/5),
    ])
def test_srgb_gamut(colorspace, n=10):
    colorio.show_srgb_gamut(colorspace, 'srgb.vtu', n=n)
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


if __name__ == '__main__':
    # colorspace_ = colorio.SrgbLinear()
    # colorspace_ = colorio.XYZ()
    # colorspace_ = colorio.XYY()
    # colorspace_ = colorio.JzAzBz()
    # colorspace_ = colorio.CIELUV()
    # colorspace_ = colorio.CIELAB()
    # colorspace_ = colorio.CAM02('UCS', 0.69, 20, 64/numpy.pi/5)
    colorspace_ = colorio.CAM16UCS(0.69, 20, 64/numpy.pi/5)
    test_srgb_gamut(colorspace_, n=50)
