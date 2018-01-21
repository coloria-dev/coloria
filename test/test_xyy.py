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
    xyy = colorio.XYY()
    out = xyy.to_xyz(xyy.from_xyz(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


def test_plot():
    colorio.XYY().srgb_gamut(n=10)
    return


def test_gamut_diagram():
    colorio.xyy.show_gamut_diagram()
    return


if __name__ == '__main__':
    test_gamut_diagram()


if __name__ == '__main__':
    test_plot()
