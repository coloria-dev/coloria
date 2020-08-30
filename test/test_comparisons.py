# import colour
import colorspacious
import numpy
from cam16_legacy import CAM16Legacy

import colorio

numpy.random.seed(0)


def test_0():
    Y_b = 20
    L_A = 64 / numpy.pi / 5
    c = 0.69  # average
    cam16 = colorio.CAM16(c, Y_b, L_A)

    xyz = numpy.zeros(3)
    J, C, _, h, M, s, Q = cam16.from_xyz100(xyz)

    assert J == 0.0
    assert C == 0.0
    assert h == 0.0
    assert M == 0.0
    assert s == 0.0
    assert Q == 0.0

    # Comparison with other schemes
    cam16_legacy = CAM16Legacy(c, Y_b, L_A)
    ref2 = cam16_legacy.from_xyz100(xyz)
    print(ref2)
    return


def test_from():
    """Compare colorio with colorspacius and colour."""
    xyz = 100 * numpy.random.rand(3)

    Y_b = 20
    whitepoint = colorio.illuminants.whitepoints_cie1931["D65"]
    L_A = 64 / numpy.pi / 5

    c = 0.69  # average
    cs2 = colorio.CIECAM02(c, Y_b, L_A)
    J, C, H, h, M, s, Q = cs2.from_xyz100(xyz)

    # compare with colorspacious
    cs1 = colorspacious.ciecam02.CIECAM02Space(
        whitepoint, Y_b, L_A, surround=colorspacious.CIECAM02Surround.AVERAGE
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
    # TODO reinstate once https://github.com/colour-science/colour/issues/467 is fixed
    # ref2 = colour.appearance.ciecam02.XYZ_to_CIECAM02(xyz, whitepoint, L_A, Y_b)
    # assert abs(ref2.J - J) < 1.0e-14 * J
    # assert abs(ref2.C - C) < 1.0e-14 * C
    # # assert abs(ref2.H - H) < 1.0e-14 * H
    # assert abs(ref2.h - h) < 1.0e-14 * h
    # assert abs(ref2.M - M) < 1.0e-14 * M
    # assert abs(ref2.s - s) < 1.0e-14 * s
    # assert abs(ref2.Q - Q) < 1.0e-14 * Q
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
    L_A = 64 / numpy.pi / 5
    c = 0.69  # average
    cam16 = colorio.CAM16(c, Y_b, L_A)

    cam16_legacy = CAM16Legacy(c, Y_b, L_A)

    perfplot.show(
        setup=setup,
        kernels=[cam16.from_xyz100, cam16_legacy.from_xyz100],
        labels=["new", "legacy"],
        n_range=1000 * numpy.arange(6),
        equality_check=False,
    )
    return


def performance_comparison_to():
    import perfplot

    Y_b = 20
    L_A = 64 / numpy.pi / 5

    c = 0.69  # average
    cam16 = colorio.CAM16(c, Y_b, L_A)

    def cio(x):
        return cam16.to_xyz100(x, "JCh")

    cam16_legacy = CAM16Legacy(c, Y_b, L_A)

    def cio_legacy(x):
        return cam16_legacy.to_xyz100(x, "JCh")

    perfplot.plot(
        setup=lambda n: numpy.random.rand(3, n),
        kernels=[cio, cio_legacy],
        n_range=100000 * numpy.arange(11),
        xlabel="Number of input samples",
    )

    # import matplotlib2tikz
    # matplotlib2tikz.save('fig.tikz')
    return


if __name__ == "__main__":
    performance_comparison_to()
    # test_0()
