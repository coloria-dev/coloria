# -*- coding: utf-8 -*-
#
import numpy
import pytest

import colorio


@pytest.mark.parametrize('vals', [
    numpy.random.rand(3),
    numpy.random.rand(3, 7),
    ])
def test_conversion(vals):
    srgb_linear = colorio.SrgbLinear()

    out = srgb_linear.to_xyz(srgb_linear.from_xyz(vals))
    assert numpy.all(abs(vals - out) < 1.0e-14)

    out = srgb_linear.to_srgb1(srgb_linear.from_srgb1(vals))
    assert numpy.all(abs(vals - out) < 1.0e-14)
    return
