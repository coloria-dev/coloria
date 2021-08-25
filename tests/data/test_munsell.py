import pytest

import colorio


def test_show():
    cs = colorio.cs.IPT
    plt = colorio.data.Munsell().plot(cs, V=5)
    plt.show()
    plt.close()


def test_show_lightness():
    cs = colorio.cs.IPT
    plt = colorio.data.Munsell().plot_lightness(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.CIELAB, 0.6928265531073651),
        (colorio.cs.CIELUV, 0.6928265531073651),
        (colorio.cs.CIELCH, 0.6928265531073651),
        (colorio.cs.CIEHCL, 0.6928265531073651),
        (colorio.cs.SRLAB2, 1.0538888582002768),
        (colorio.cs.RLAB, 2.7274797600841056),
        (colorio.cs.CAM02SCD, 4.342489913211733),
        (colorio.cs.CAM02UCS, 4.342489913211733),
        (colorio.cs.CAM02LCD, 4.342489913211733),
        (colorio.cs.CAM16UCS, 4.704019435432123),
        (colorio.cs.JzAzBz, 8.90523847167804),
        (colorio.cs.IPT, 9.47504412245566),
        (colorio.cs.OsaUcs, 10.706056005621226),
        (colorio.cs.OKLAB, 11.192466987558433),
        (colorio.cs.PROLAB, 11.81141591894177),
        (colorio.cs.XYY100, 31.251745570029676),
        # (colorio.cs.ICtCp, nan),
    ],
)
def test_stress(cs_class, ref):
    print(cs_class)
    res = colorio.data.Munsell().stress_lightness(cs_class)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    # test_stress()
    # test_show()
    test_show_lightness()
