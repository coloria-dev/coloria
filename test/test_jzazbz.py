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
    jzazbz = colorio.JzAzBz()
    out = jzazbz.to_xyz100(jzazbz.from_xyz100(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-13)
    return
