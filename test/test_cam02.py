# -*- coding: utf-8 -*-
#
import numpy
import pytest

import colorio

numpy.random.seed(0)


@pytest.mark.parametrize('variant', [
    'LCD', 'SCD', 'UCS',
    ])
@pytest.mark.parametrize('xyz', [
    numpy.random.rand(3),
    numpy.random.rand(3, 7),
    ])
def test_conversion(variant, xyz):
    # test with srgb conditions
    L_A = 64 / numpy.pi / 5
    cam02 = colorio.CAM02(variant, 0.69, 20, L_A)
    out = cam02.to_xyz(cam02.from_xyz(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return
