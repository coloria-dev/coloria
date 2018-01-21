# -*- coding: utf-8 -*-
#
import numpy
import pytest

import colorio

numpy.random.seed(0)


@pytest.mark.parametrize('xyz', [
    numpy.random.rand(3),
    numpy.random.rand(3, 7),
    ])
def test_conversion(xyz):
    # test with srgb conditions
    L_A = 64 / numpy.pi / 5
    cam02_ucs = colorio.CAM02_UCS(0.69, 20, L_A)
    out = cam02_ucs.to_xyz(cam02_ucs.from_xyz(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


def test_srgb_gamut(n=10):
    L_A = 64 / numpy.pi / 5
    cam02_ucs = colorio.CAM02_UCS(0.69, 20, L_A)
    cam02_ucs.srgb_gamut(n=n)
    return


if __name__ == '__main__':
    test_srgb_gamut(n=50)
