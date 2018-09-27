# -*- coding: utf-8 -*-
#
import numpy
import pytest

import colorio

numpy.random.seed(0)


@pytest.mark.parametrize(
    "xyz", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion(xyz):
    # test with srgb conditions
    L_A = 64 / numpy.pi / 5
    ciecam02 = colorio.CIECAM02(0.69, 20, L_A)
    J, C, H, h, M, s, Q = ciecam02.from_xyz100(xyz)

    out = ciecam02.to_xyz100(numpy.array([J, C, H]), "JCH")
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = ciecam02.to_xyz100(numpy.array([Q, M, h]), "QMh")
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))

    out = ciecam02.to_xyz100(numpy.array([J, s, h]), "Jsh")
    assert numpy.all(abs(xyz - out) < 1.0e-13 * abs(xyz))
    return


def test_breakdown():
    bad_xyz = [8.71292997, 2.02183974, 83.26198455]
    L_A = 64 / numpy.pi / 5
    ciecam02 = colorio.CIECAM02(0.69, 20, L_A)
    with pytest.raises(colorio.ciecam02.NegativeAError):
        ciecam02.from_xyz100(bad_xyz)
    return


@pytest.mark.parametrize("variant", ["LCD", "SCD", "UCS"])
@pytest.mark.parametrize(
    "xyz", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion_variants(variant, xyz):
    # test with srgb conditions
    L_A = 64 / numpy.pi / 5
    cam02 = colorio.CAM02(variant, 0.69, 20, L_A)
    out = cam02.to_xyz100(cam02.from_xyz100(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-14)
    return


@pytest.mark.parametrize("xyz", [numpy.zeros(3), numpy.zeros((3, 4, 5))])
def test_zero(xyz):
    L_A = 64 / numpy.pi / 5
    cs = colorio.CIECAM02(0.69, 20, L_A)
    J, C, H, h, M, s, Q = cs.from_xyz100(xyz)

    assert numpy.all(J == 0.0)
    assert numpy.all(C == 0.0)
    assert numpy.all(h == 0.0)
    assert numpy.all(M == 0.0)
    assert numpy.all(s == 0.0)
    assert numpy.all(Q == 0.0)

    out = cs.to_xyz100(numpy.array([J, C, H]), "JCH")
    assert numpy.all(abs(out) < 1.0e-13)

    out = cs.to_xyz100(numpy.array([Q, M, h]), "QMh")
    assert numpy.all(abs(out) < 1.0e-13)

    out = cs.to_xyz100(numpy.array([J, s, h]), "Jsh")
    assert numpy.all(abs(out) < 1.0e-13)
    return


def test_gold():
    # See
    # https://github.com/njsmith/colorspacious/blob/master/colorspacious/gold_values.py
    xyz = [19.31, 23.93, 10.14]
    ciecam02 = colorio.CIECAM02(0.69, 18, 200, whitepoint=[98.88, 90, 32.03])
    values = ciecam02.from_xyz100(xyz)
    reference_values = numpy.array(
        [
            # J, C, H, h, M, s, Q
            48.0314,
            38.7789,
            240.8885,
            191.0452,
            38.7789,
            46.0177,
            183.1240,
        ]
    )
    assert numpy.all(abs(values - reference_values) < 1.0e-6 * reference_values)

    # different L_A
    xyz = [19.31, 23.93, 10.14]
    ciecam02 = colorio.CIECAM02(0.69, 18, 20, whitepoint=[98.88, 90, 32.03])
    values = ciecam02.from_xyz100(xyz)
    reference_values = numpy.array(
        [
            # J, C, H, h, M, s, Q
            47.6856,
            36.0527,
            232.6630,
            185.3445,
            29.7580,
            51.1275,
            113.8401,
        ]
    )
    assert numpy.all(abs(values - reference_values) < 1.0e-5 * reference_values)

    # gold values from Mark Fairchild's spreadsheet at
    #   http://rit-mcsl.org/fairchild//files/AppModEx.xls
    xyz = [19.01, 20.00, 21.78]
    ciecam02 = colorio.CIECAM02(
        0.69, 20, 318.30988618379, whitepoint=[95.05, 100.0, 108.88]
    )
    values = ciecam02.from_xyz100(xyz)
    reference_values = numpy.array(
        [
            # J, C, H, h, M, s, Q
            4.17310911e01,
            1.04707861e-01,
            2.78060724e02,
            2.19048423e02,
            1.08842280e-01,
            2.36030659e00,
            1.95371311e02,
        ]
    )
    assert numpy.all(abs(values - reference_values) < 1.0e-8 * reference_values)

    xyz = [57.06, 43.06, 31.96]
    ciecam02 = colorio.CIECAM02(
        0.69, 20, 31.830988618379, whitepoint=[95.05, 100.0, 108.88]
    )
    values = ciecam02.from_xyz100(xyz)
    reference_values = numpy.array(
        [
            # J, C, H, h, M, s, Q
            65.95523,
            48.57050,
            399.38837,
            19.55739,
            41.67327,
            52.24549,
            152.67220,
        ]
    )
    assert numpy.all(abs(values - reference_values) < 1.0e-6 * reference_values)
    return


# @pytest.mark.parametrize('variant, xyz100, ref', [
#     # From
# <https://github.com/njsmith/colorspacious/blob/master/colorspacious/gold_values.py>.
#     ('UCS', [50, 20, 10], [62.96296296, 16.22742674, 2.86133316]),
#     ('UCS', [10, 60, 100], [15.88785047, -6.56546789, 37.23461867]),
#     ('LCD', [50, 20, 10], [81.77008177, 18.72061994, 3.30095039]),
#     ('LCD', [10, 60, 100], [20.63357204, -9.04659289, 51.30577777]),
#     ('SCD', [50, 20, 10], [50.77658303, 14.80756375, 2.61097301]),
#     ('SCD', [10, 60, 100], [12.81278263, -5.5311588, 31.36876036]),
#     ])
# def test_reference(variant, xyz100, ref):
#     cs = colorio.CAM02(
#         variant, whitepoint=numpy.array([96.422, 100, 82.521])/100
#         )
#     xyz = numpy.array(xyz100) / 100
#     assert numpy.all(
#         abs(cs.from_xyz100(xyz) - ref) < 1.0e-6 * abs(numpy.array(ref))
#         )
#     return


if __name__ == "__main__":
    test_gold()
