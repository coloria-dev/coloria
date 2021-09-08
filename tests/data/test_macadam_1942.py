import pytest

import colorio


def test_show():
    cs = colorio.cs.XYY1
    plt = colorio.data.MacAdam1942(Y=50.0).plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.CAM02LCD, 26.984474619627647),
        (colorio.cs.CAM02UCS, 31.831373461580487),
        (colorio.cs.CAM16UCS, 34.00812991052416),
        (colorio.cs.CAM02SCD, 35.988912532429346),
        (colorio.cs.CIELAB, 44.95432696577654),
        (colorio.cs.CIELCH, 62.53845922365119),
        (colorio.cs.CIELUV, 36.50102760327483),
        (colorio.cs.CIEHCL, 66.31832189088743),
        (colorio.cs.ICtCp, 77.23169033651011),
        (colorio.cs.IPT, 38.4833636095126),
        (colorio.cs.JzAzBz, 33.27570881193104),
        (colorio.cs.OKLAB, 37.67725482006628),
        (colorio.cs.OsaUcs, 42.43966467973784),
        (colorio.cs.PROLAB, 35.105986390925786),
        (colorio.cs.RLAB, 45.99760605766801),
        (colorio.cs.SRLAB2, 39.19981054654803),
        (colorio.cs.XYY100, 53.750602251102265),
        (colorio.cs.XYZ100, 79.56476212916836),
    ],
)
def test_stress(cs_class, ref):
    print(cs_class)
    res = colorio.data.MacAdam1942(50.0).stress(cs_class)
    print(res)
    assert abs(res - ref) < 1.0e-13 * ref


if __name__ == "__main__":
    test_show()
    # test_stress()
