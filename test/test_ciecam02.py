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
    ciecam02 = colorio.CIECAM02(0.69, 20, L_A)
    J, C, H, h, M, s, Q = ciecam02.from_xyz(xyz)

    out = ciecam02.to_xyz(numpy.array([J, C, H]), 'JCH')
    assert numpy.all(abs(xyz - out) < 1.0e-14)

    out = ciecam02.to_xyz(numpy.array([Q, M, h]), 'QMh')
    assert numpy.all(abs(xyz - out) < 1.0e-14)

    out = ciecam02.to_xyz(numpy.array([J, s, h]), 'Jsh')
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


def test_breakdown():
    bad_xyz = [8.71292997, 2.02183974, 83.26198455]
    L_A = 64 / numpy.pi / 5
    ciecam02 = colorio.CIECAM02(0.69, 20, L_A)
    with pytest.raises(colorio.ciecam02.NegativeAError):
        ciecam02.from_xyz(bad_xyz)
    return


@pytest.mark.parametrize('variant', [
    'LCD', 'SCD', 'UCS',
    ])
@pytest.mark.parametrize('xyz', [
    numpy.random.rand(3),
    numpy.random.rand(3, 7),
    ])
def test_conversion_variants(variant, xyz):
    # test with srgb conditions
    L_A = 64 / numpy.pi / 5
    cam02 = colorio.CAM02(variant, 0.69, 20, L_A)
    out = cam02.to_xyz(cam02.from_xyz(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


if __name__ == '__main__':
    test_breakdown()
