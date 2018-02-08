# -*- coding: utf-8 -*-
#
import numpy

import colour
import colorio
import colorspacious

numpy.random.seed(0)


def test_0():
    Y_b = 20
    whitepoint = colorio.illuminants.whitepoints_cie1931['D65']
    L_A = 64 / numpy.pi / 5

    c = 0.69  # average
    cs2 = colorio.CIECAM02(c, Y_b, L_A)

    xyz = numpy.zeros(3)
    J, C, _, h, M, s, Q = cs2.from_xyz100(xyz)

    assert J == 0.0
    assert C == 0.0
    assert h == 0.0
    assert M == 0.0
    assert s == 0.0
    assert Q == 0.0

    # Comparison with other schemes
    cs1 = colorspacious.ciecam02.CIECAM02Space(
        whitepoint, Y_b, L_A,
        surround=colorspacious.CIECAM02Surround.AVERAGE
        )
    ref1 = cs1.XYZ100_to_CIECAM02(xyz)
    print(ref1)

    ref2 = colour.appearance.ciecam02.XYZ_to_CIECAM02(
        xyz, whitepoint, L_A, Y_b
        )
    print(ref2)
    return


def test_from():
    '''Compare colorio with colorspacius and colour.
    '''
    xyz = 100 * numpy.random.rand(3)

    Y_b = 20
    whitepoint = colorio.illuminants.whitepoints_cie1931['D65']
    L_A = 64 / numpy.pi / 5

    c = 0.69  # average
    cs2 = colorio.CIECAM02(c, Y_b, L_A)
    J, C, H, h, M, s, Q = cs2.from_xyz100(xyz)

    # compare with colorspacious
    cs1 = colorspacious.ciecam02.CIECAM02Space(
        whitepoint, Y_b, L_A,
        surround=colorspacious.CIECAM02Surround.AVERAGE
        )
    ref1 = cs1.XYZ100_to_CIECAM02(xyz)
    assert abs(ref1.J - J) < 1.0e-14 * J
    assert abs(ref1.C - C) < 1.0e-14 * C
    assert abs(ref1.H - H) < 1.0e-14 * H
    assert abs(ref1.h - h) < 1.0e-14 * h
    assert abs(ref1.M - M) < 1.0e-14 * M
    assert abs(ref1.s - s) < 1.0e-14 * s
    assert abs(ref1.Q - Q) < 1.0e-14 * Q

    # compare with color
    ref2 = colour.appearance.ciecam02.XYZ_to_CIECAM02(
        xyz, whitepoint, L_A, Y_b
        )
    assert abs(ref2.J - J) < 1.0e-14 * J
    assert abs(ref2.C - C) < 1.0e-14 * C
    assert abs(ref2.H - H) < 1.0e-14 * H
    assert abs(ref2.h - h) < 1.0e-14 * h
    assert abs(ref2.M - M) < 1.0e-14 * M
    assert abs(ref2.s - s) < 1.0e-14 * s
    assert abs(ref2.Q - Q) < 1.0e-14 * Q
    return


def performance_comparison_from():
    import perfplot

    def setup(n):
        out = numpy.empty((3, n))
        rgb = numpy.random.rand(3)
        for k in range(3):
            out[k] = rgb[k]
        return out

    Y_b = 20
    whitepoint = colorio.illuminants.whitepoints_cie1931['D65']
    L_A = 64 / numpy.pi / 5
    cs1 = colorspacious.ciecam02.CIECAM02Space(
        whitepoint, Y_b, L_A,
        surround=colorspacious.CIECAM02Surround.AVERAGE
        )

    c = 0.69  # average
    cs2 = colorio.CIECAM02(c, Y_b, L_A)

    perfplot.show(
        setup=setup,
        kernels=[
            lambda x: cs1.XYZ100_to_CIECAM02(x.T),
            cs2.from_xyz100,
            lambda x: colour.appearance.ciecam02.XYZ_to_CIECAM02(
                x.T, whitepoint, L_A, Y_b
                )
            ],
        labels=['colorspacious', 'colorio', 'colour'],
        n_range=10000 * numpy.arange(6),
        equality_check=False
        )
    return


def performance_comparison_to():
    import perfplot

    Y_b = 20
    whitepoint = colorio.illuminants.whitepoints_cie1931['D65']
    L_A = 64 / numpy.pi / 5

    c = 0.69  # average
    cs2 = colorio.CIECAM02(c, Y_b, L_A)

    cs1 = colorspacious.ciecam02.CIECAM02Space(
        whitepoint, Y_b, L_A,
        surround=colorspacious.CIECAM02Surround.AVERAGE
        )

    def setup(n):
        rgb = numpy.random.rand(3)
        out = numpy.empty((3, n))
        for k in range(3):
            out[k] = rgb[k]
        return out

    def csp(x):
        return cs1.CIECAM02_to_XYZ100(J=x[0], C=x[1], h=x[2])

    def cio(x):
        return cs2.to_xyz100(x, 'JCh')

    def clr(x):
        J, C, h = x
        spec = colour.appearance.ciecam02.CIECAM02_Specification(J=J, C=C, h=h)
        return colour.appearance.ciecam02.CIECAM02_to_XYZ(
                spec, whitepoint, L_A, Y_b
                )

    perfplot.show(
        setup=setup,
        kernels=[
            # csp,
            cio, clr
            ],
        n_range=100000 * numpy.arange(11),
        equality_check=False,
        xlabel='Number of input samples'
        )

    # import matplotlib2tikz
    # matplotlib2tikz.save('fig.tikz')
    return


if __name__ == '__main__':
    performance_comparison_to()
    # test_0()
