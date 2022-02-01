import pytest

import colorio


def test_show():
    cs = colorio.cs.CIELAB
    plt = colorio.data.RitDupont().plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.CAM16UCS, 20.56944422527809),
        (colorio.cs.CAM02UCS, 21.272923860475228),
        (colorio.cs.CAM02SCD, 21.886793024730252),
        (colorio.cs.CAM02LCD, 26.104986244230254),
        (colorio.cs.JzAzBz, 28.29909799913785),
        (colorio.cs.IPT, 29.447144074050673),
        (colorio.cs.PROLAB, 30.168608588989056),
        (colorio.cs.OsaUcs, 31.486420519127105),
        (colorio.cs.OKLAB, 31.81741455159825),
        (colorio.cs.SRLAB2, 32.798157162944854),
        (colorio.cs.CIELAB, 33.41593673343391),
        (colorio.cs.RLAB, 33.986144075747816),
        (colorio.cs.CIELUV, 36.529835101856975),
        (colorio.cs.ICtCp, 50.812320953414655),
        (colorio.cs.XYZ100, 61.03890586028899),
        (colorio.cs.XYY100, 72.02049789265273),
        (colorio.cs.CIEHCL, 96.11444826150097),
        (colorio.cs.CIELCH, 96.17992567400712),
    ],
)
def test_stress(cs_class, ref):
    res = colorio.data.RitDupont().stress(cs_class)
    print(cs_class)
    print(res)
    assert abs(res - ref) < 1.0e-13 * ref


@pytest.mark.parametrize(
    "fun,ref",
    [
        (colorio.diff.ciede2000, 19.46979726886805),
        (colorio.diff.cie94, 20.29959973965712),
        (colorio.diff.cmc, 33.208917491171825),
        (colorio.diff.cie76, 33.41593673343391),
    ],
)
def test_stress_diff(fun, ref):
    print(fun)
    data = colorio.data.RitDupont()
    res = data.stress_lab_diff(fun)
    print(res)
    assert abs(res - ref) < 1.0e-12 * ref
