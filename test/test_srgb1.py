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
    out = colorio.srgb1.to_srgb_linear(colorio.srgb1.from_srgb_linear(vals))
    print(vals)
    print(out)
    assert numpy.all(abs(vals - out) < 1.0e-14)
    return


def test_gamut():
    colorio.srgb1.srgb_gamut(n=10)
    return


if __name__ == '__main__':
    test_gamut()
