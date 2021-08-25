import pytest

import colorio


def test_show():
    cs = colorio.cs.CAM02UCS
    plt = colorio.data.Witt().plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.CAM02SCD, 28.577269461047187),
        (colorio.cs.CAM02UCS, 30.4584731112042),
        (colorio.cs.CAM16UCS, 30.991751511774744),
        (colorio.cs.CAM02LCD, 39.16831128020275),
        (colorio.cs.OKLAB, 45.2276225067436),
        (colorio.cs.IPT, 46.483167491325595),
        (colorio.cs.OsaUcs, 47.3978125741753),
        (colorio.cs.JzAzBz, 48.53002773753702),
        (colorio.cs.RLAB, 49.87849727183007),
        (colorio.cs.SRLAB2, 51.549289192015046),
        (colorio.cs.CIELAB, 51.70889935711488),
        (colorio.cs.PROLAB, 51.8504406525641),
        (colorio.cs.XYZ100, 52.33740794765402),
        (colorio.cs.CIELUV, 53.200886880239736),
        (colorio.cs.ICtCp, 60.43870870332812),
        (colorio.cs.XYY100, 77.12352248970204),
        (colorio.cs.CIELCH, 87.43172873664071),
        (colorio.cs.CIEHCL, 87.57481132806802),
    ],
)
def test_stress(cs_class, ref):
    res = colorio.data.Witt().stress(cs_class)
    print(cs_class)
    print(res)
    assert abs(res - ref) < 1.0e-11 * ref


@pytest.mark.parametrize(
    "fun,ref",
    [
        (colorio.diff.ciede2000, 30.218234357300684),
        (colorio.diff.cie94, 31.704895635913026),
        (colorio.diff.cmc, 42.1795697284095),
        (colorio.diff.cie76, 51.70889935711488),
    ],
)
def test_stress_diff(fun, ref):
    print(fun)
    data = colorio.data.Witt()
    res = data.stress_lab_diff(fun)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref
