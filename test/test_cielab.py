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
    print(xyz)
    out = colorio.cielab.to_xyz(colorio.cielab.from_xyz(xyz))
    print(out)
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


# def test_luminance_level():
#     colorio.cielab.show_luminance_level(50)
#     return


if __name__ == '__main__':
    colorio.cielab.show_luminance_level(50)
    # test_luminance_level()
