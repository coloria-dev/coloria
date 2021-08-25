import pytest

import colorio


def test_show():
    cs = colorio.cs.JzAzBz
    plt = colorio.data.FairchildChen("SL2").plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref1,ref2",
    [
        (colorio.cs.JzAzBz, 4.079328782823399, 10.428911813807975),
        (colorio.cs.SRLAB2, 4.669640386357883, 11.272597307546123),
        (colorio.cs.CIELUV, 4.672891455041493, 11.266510241793169),
        (colorio.cs.CIELAB, 4.672891455041493, 11.266510241793169),
        (colorio.cs.CIELCH, 4.672891455041493, 11.266510241793169),
        (colorio.cs.CIEHCL, 4.672891455041493, 11.266510241793169),
        (colorio.cs.IPT, 5.1650760027712534, 9.428338645588513),
        (colorio.cs.RLAB, 5.277504252399791, 9.44336537734861),
        (colorio.cs.OKLAB, 5.646655086765461, 13.12276913806233),
        (colorio.cs.OsaUcs, 5.827386521775151, 13.443081202978794),
        (colorio.cs.CAM16UCS, 18.64772022812977, 29.729084277208024),
        (colorio.cs.CAM02SCD, 18.649400893925677, 29.75664952154394),
        (colorio.cs.CAM02UCS, 18.649400893925677, 29.75664952154394),
        (colorio.cs.CAM02LCD, 18.649400893925677, 29.75664952154394),
        (colorio.cs.ICtCp, 25.897107754911293, 39.360406514124904),
        (colorio.cs.PROLAB, 29.984387363692317, 43.83803651177642),
        (colorio.cs.XYY100, 35.9945312574285, 41.68242289798867),
    ],
)
def test_stress(cs_class, ref1, ref2):
    print(cs_class)
    data = colorio.data.FairchildChen("SL1")
    res = data.stress(cs_class)
    print(res)
    assert abs(res - ref1) < 1.0e-14 * ref1

    data = colorio.data.FairchildChen("SL2")
    res = data.stress(cs_class)
    print(res)
    assert abs(res - ref2) < 1.0e-14 * ref2


if __name__ == "__main__":
    test_show()
