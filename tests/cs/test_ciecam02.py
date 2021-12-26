import numpy as np
import pytest

import colorio

rng = np.random.default_rng(0)


@pytest.mark.parametrize(
    "xyz", [rng.random(3), rng.random((3, 7)), rng.random((3, 4, 5))]
)
def test_conversion(xyz, tol=1.0e-12):
    # test with srgb conditions
    ciecam02 = colorio.cs.CIECAM02(0.69, 20, 20)
    J, C, H, h, M, s, Q = ciecam02.from_xyz100(xyz)

    out = ciecam02.to_xyz100(np.array([J, C, H]), "JCH")
    assert np.all(np.abs(xyz - out) < tol * np.abs(xyz))

    out = ciecam02.to_xyz100(np.array([Q, M, h]), "QMh")
    assert np.all(np.abs(xyz - out) < tol * np.abs(xyz))

    out = ciecam02.to_xyz100(np.array([J, s, h]), "Jsh")
    assert np.all(np.abs(xyz - out) < tol * np.abs(xyz))


def test_breakdown():
    bad_xyz = [8.71292997, 2.02183974, 83.26198455]
    ciecam02 = colorio.cs.CIECAM02(0.69, 20, 20)
    with pytest.raises(colorio.ColorioError):
        ciecam02.from_xyz100(bad_xyz)


@pytest.mark.parametrize("xyz", [np.zeros(3), np.zeros((3, 4, 5))])
def test_zero(xyz):
    cs = colorio.cs.CIECAM02(0.69, 20, 20)
    J, C, H, h, M, s, Q = cs.from_xyz100(xyz)

    assert np.all(J == 0.0)
    assert np.all(C == 0.0)
    assert np.all(h == 0.0)
    assert np.all(M == 0.0)
    assert np.all(s == 0.0)
    assert np.all(Q == 0.0)

    out = cs.to_xyz100(np.array([J, C, H]), "JCH")
    assert np.all(abs(out) < 1.0e-13)

    out = cs.to_xyz100(np.array([Q, M, h]), "QMh")
    assert np.all(abs(out) < 1.0e-13)

    out = cs.to_xyz100(np.array([J, s, h]), "Jsh")
    assert np.all(abs(out) < 1.0e-13)


def test_gold():
    # See
    # https://github.com/njsmith/colorspacious/blob/master/colorspacious/gold_values.py
    xyz = [19.31, 23.93, 10.14]
    ciecam02 = colorio.cs.CIECAM02(0.69, 18, 200, whitepoint=[98.88, 90, 32.03])
    values = ciecam02.from_xyz100(xyz)
    # J, C, H, h, M, s, Q
    reference_values = [
        48.0314,
        38.7789,
        240.8884,
        191.0452,
        38.7789,
        46.0177,
        183.1240,
    ]
    assert np.all(values.round(4) == reference_values)

    # different L_A
    xyz = [19.31, 23.93, 10.14]
    ciecam02 = colorio.cs.CIECAM02(0.69, 18, 20, whitepoint=[98.88, 90, 32.03])
    values = ciecam02.from_xyz100(xyz)
    # J, C, H, h, M, s, Q
    reference_values = [
        47.6856,
        36.0527,
        232.6630,
        185.3445,
        29.7580,
        51.1275,
        113.8401,
    ]
    assert np.all(values.round(4) == reference_values)

    # gold values from Mark Fairchild's spreadsheet at
    # http://rit-mcsl.org/fairchild//files/AppModEx.xls
    xyz = [19.01, 20.00, 21.78]
    ciecam02 = colorio.cs.CIECAM02(
        0.69, 20, 318.30988618379, whitepoint=[95.05, 100.0, 108.88]
    )
    values = ciecam02.from_xyz100(xyz)
    reference_values = np.array(
        # J, C, H, h, M, s, Q
        [
            4.17310911e01,
            1.04707861e-01,
            2.78060724e02,
            2.19048423e02,
            1.08842280e-01,
            2.36030659e00,
            1.95371311e02,
        ]
    )
    assert np.all(np.abs(values - reference_values) < 1.0e-8 * reference_values)

    xyz = [57.06, 43.06, 31.96]
    ciecam02 = colorio.cs.CIECAM02(
        0.69, 20, 31.830988618379, whitepoint=[95.05, 100.0, 108.88]
    )
    values = ciecam02.from_xyz100(xyz)
    # J, C, H, h, M, s, Q
    reference_values = [
        65.95523,
        48.57048,
        399.38844,
        19.55738,
        41.67325,
        52.24548,
        152.67220,
    ]
    print(values.round(5))
    assert np.all(values.round(5) == reference_values)


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
#         variant, whitepoint=np.array([96.422, 100, 82.521])/100
#         )
#     xyz = np.array(xyz100) / 100
#     assert np.all(
#         abs(cs.from_xyz100(xyz) - ref) < 1.0e-6 * abs(np.array(ref))
#         )
