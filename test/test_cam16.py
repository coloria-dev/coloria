# -*- coding: utf-8 -*-
#
import numpy
import pytest

import colorio

numpy.random.seed(0)


@pytest.mark.parametrize('xyz', [
    100 * numpy.random.rand(3),
    100 * numpy.random.rand(3, 7),
    100 * numpy.random.rand(3, 4, 5),
    ])
def test_conversion(xyz):
    # test with srgb conditions
    L_A = 64 / numpy.pi / 5
    cam16 = colorio.CAM16(0.69, 20, L_A)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    out = cam16.to_xyz100(numpy.array([J, C, H]), 'JCH')
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = cam16.to_xyz100(numpy.array([Q, M, h]), 'QMh')
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = cam16.to_xyz100(numpy.array([J, s, h]), 'Jsh')
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))
    return


def test_zero():
    xyz = numpy.zeros(3)

    L_A = 64 / numpy.pi / 5
    cam16 = colorio.CAM16(0.69, 20, L_A)
    J, C, H, h, M, s, Q = cam16.from_xyz100(xyz)

    assert J == 0.0
    assert C == 0.0
    assert h == 0.0
    assert M == 0.0
    assert s == 0.0
    assert Q == 0.0

    out = cam16.to_xyz100(numpy.array([J, C, H]), 'JCH')
    # assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))
    print(out)
    exit(1)
    return


@pytest.mark.parametrize('xyz', [
    numpy.random.rand(3),
    numpy.random.rand(3, 7),
    numpy.random.rand(3, 4, 5),
    ])
def test_conversion_variants(xyz):
    # test with srgb conditions
    L_A = 64 / numpy.pi / 5
    cam16 = colorio.CAM16UCS(0.69, 20, L_A)
    out = cam16.to_xyz100(cam16.from_xyz100(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


if __name__ == '__main__':
    test_conversion(numpy.random.rand(3))
