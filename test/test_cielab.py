# -*- coding: utf-8 -*-
#
import numpy
import pytest

import colorio


@pytest.mark.parametrize('xyz', [
    numpy.random.rand(3),
    numpy.random.rand(3, 7),
    ])
def test_conversion(xyz):
    cielab = colorio.CIELAB()
    out = cielab.to_xyz(cielab.from_xyz(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


def test_srgb_gamut():
    colorio.CIELAB().srgb_gamut(n=10)
    return


if __name__ == '__main__':
    test_srgb_gamut()
