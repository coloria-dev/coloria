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
    cam16 = colorio.CAM16(0.69, 20, L_A)
    J, C, H, h, M, s, Q = cam16.from_xyz(xyz)

    out = cam16.to_xyz(numpy.array([J, C, H]), 'JCH')
    assert numpy.all(abs(xyz - out) < 1.0e-14)

    out = cam16.to_xyz(numpy.array([Q, M, h]), 'QMh')
    assert numpy.all(abs(xyz - out) < 1.0e-14)

    out = cam16.to_xyz(numpy.array([J, s, h]), 'Jsh')
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


@pytest.mark.parametrize('xyz', [
    numpy.random.rand(3),
    numpy.random.rand(3, 7),
    ])
def test_conversion_variants(xyz):
    # test with srgb conditions
    L_A = 64 / numpy.pi / 5
    cam16 = colorio.CAM16UCS(0.69, 20, L_A)
    out = cam16.to_xyz(cam16.from_xyz(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


if __name__ == '__main__':
    test_conversion(numpy.random.rand(3))
