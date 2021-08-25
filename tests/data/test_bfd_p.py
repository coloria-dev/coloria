import pytest

import colorio

# def test_show():
#     # cs = colorio.cs.CIELAB()
#     # cs = colorio.cs.CIEHCL()
#     # cs = colorio.cs.CIELCH()
#     # cs = colorio.cs.OsaUcs()
#     # cs = colorio.cs.IPT()
#     # cs = colorio.cs.OKLAB()
#     cs = colorio.cs.CAM02("UCS", 0.59, 20, 20)
#     # cs = colorio.cs.CAM16UCS(0.59, 20, 20)
#     # cs = colorio.cs.JzAzBz()
#     # cs = colorio.cs.XYY(1)
#     plt = colorio.data.BfdP().plot(cs)
#     plt.show()
#     plt.close()


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.CAM16UCS, 30.912784895833845),
        (colorio.cs.CAM02UCS, 30.96274433209427),
        (colorio.cs.CAM02SCD, 32.74387216723911),
        (colorio.cs.CAM02LCD, 33.69705929801087),
        (colorio.cs.JzAzBz, 37.706002260925544),
        (colorio.cs.OsaUcs, 40.16714621119245),
        (colorio.cs.IPT, 41.08905311640498),
        (colorio.cs.RLAB, 41.39045505197253),
        (colorio.cs.CIELAB, 42.46259476788477),
        (colorio.cs.SRLAB2, 42.68223438622064),
        (colorio.cs.CIELUV, 43.41087692107525),
        (colorio.cs.PROLAB, 44.40855824578013),
        (colorio.cs.OKLAB, 47.77647147881965),
        (colorio.cs.ICtCp, 55.60794680759048),
        (colorio.cs.XYZ100, 69.58673571004842),
        (colorio.cs.XYY100, 81.94385663616511),
        (colorio.cs.CIELCH, 94.1656158806),
        (colorio.cs.CIEHCL, 94.36424820618922),
    ],
)
def test_stress(cs_class, ref):
    print(cs_class)
    data = colorio.data.BfdP()
    res = data.stress(cs_class)
    print(res)
    assert abs(res - ref) < 1.0e-13 * ref


# TODO test relative
# ref = 44.17345214550509
# res = data.stress(cs, variant="relative")
# print(res)
# assert abs(res - ref) < 1.0e-14 * ref


@pytest.mark.parametrize(
    "fun,ref",
    [
        (colorio.diff.ciede2000, 29.55418997072668),
        (colorio.diff.cmc, 33.18432695639392),
        (colorio.diff.cie94, 33.70433399277904),
        (colorio.diff.cie76, 42.46259476788477),
    ],
)
def test_stress_diff(fun, ref):
    print(fun)
    data = colorio.data.BfdP()
    res = data.stress_lab_diff(fun)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref
