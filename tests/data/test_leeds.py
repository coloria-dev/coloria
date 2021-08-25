import pytest

import colorio


def test_show():
    cs = colorio.cs.CAM02UCS
    plt = colorio.data.Leeds().plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.CAM16UCS, 24.465011926041896),
        (colorio.cs.CAM02UCS, 24.594471490362963),
        (colorio.cs.CAM02SCD, 25.282643714845925),
        (colorio.cs.CAM02LCD, 30.90794447931537),
        (colorio.cs.PROLAB, 39.353476615753685),
        (colorio.cs.IPT, 39.85774625307187),
        (colorio.cs.CIELAB, 40.093223056788965),
        (colorio.cs.JzAzBz, 40.12503467515118),
        (colorio.cs.SRLAB2, 40.621076323588824),
        (colorio.cs.RLAB, 41.417149798464884),
        (colorio.cs.OsaUcs, 41.662922178924674),
        (colorio.cs.ICtCp, 43.16348979071228),
        (colorio.cs.OKLAB, 45.04796675693954),
        (colorio.cs.CIELUV, 48.29473655673663),
        (colorio.cs.XYZ100, 72.08162100159821),
        (colorio.cs.CIELCH, 78.39744306090111),
        (colorio.cs.XYY100, 79.988486428314),
        (colorio.cs.CIEHCL, 92.25243272521136),
    ],
)
def test_stress(cs_class, ref):
    print(cs_class)
    res = colorio.data.Leeds().stress(cs_class)
    print(res)
    assert abs(res - ref) < 1.0e-11 * ref


@pytest.mark.parametrize(
    "fun,ref",
    [
        (colorio.diff.ciede2000, 19.246948849846802),
        (colorio.diff.cie94, 30.494358035554026),
        (colorio.diff.cmc, 35.51506058857547),
        (colorio.diff.cie76, 40.093223056788965),
    ],
)
def test_stress_diff(fun, ref):
    print(fun)
    data = colorio.data.Leeds()
    res = data.stress_lab_diff(fun)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref
