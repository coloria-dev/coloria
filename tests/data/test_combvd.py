import pytest

import colorio


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.CAM02UCS, 29.083867343079344),
        (colorio.cs.CAM16UCS, 29.086846385696163),
        (colorio.cs.CAM02SCD, 29.995677764707914),
        (colorio.cs.CAM02LCD, 33.72500663283578),
        (colorio.cs.JzAzBz, 39.94807700406101),
        (colorio.cs.OsaUcs, 41.33380920934038),
        (colorio.cs.IPT, 41.35213361292286),
        (colorio.cs.RLAB, 43.16874979650741),
        (colorio.cs.SRLAB2, 43.96192020831966),
        (colorio.cs.CIELAB, 43.92861966265739),
        (colorio.cs.PROLAB, 44.3768306135181),
        (colorio.cs.CIELUV, 46.11433739442207),
        (colorio.cs.OKLAB, 46.12417653524034),
        (colorio.cs.ICtCp, 54.79662890386214),
        (colorio.cs.XYZ100, 67.6846121905887),
        (colorio.cs.XYY100, 80.60156328817388),
        (colorio.cs.CIELCH, 94.04196618360199),
        (colorio.cs.CIEHCL, 94.2968214535206),
    ],
)
def test_stress(cs_class, ref):
    print(cs_class)
    data = colorio.data.COMBVD()
    # ref = 43.92861966265739
    res = data.stress(cs_class)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


def test_stress_relative():
    data = colorio.data.COMBVD()
    cs = colorio.cs.CIELAB

    ref = 44.370141147103546
    res = data.stress(cs, variant="relative")
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


@pytest.mark.parametrize(
    "fun,ref",
    [
        (colorio.diff.ciede2000, 27.488683366608225),
        (colorio.diff.cie94, 31.9307958569976),
        (colorio.diff.cmc, 35.57342155223579),
        (colorio.diff.cie76, 43.92861966265739),
    ],
)
def test_stress_diff(fun, ref):
    print(fun)
    data = colorio.data.COMBVD()
    res = data.stress_lab_diff(fun)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref
