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
    out = colorio.xyy.to_xyz(colorio.xyy.from_xyz(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


def test_plot():
    colorio.xyy.srgb_gamut()
    return


if __name__ == '__main__':
    test_plot()
