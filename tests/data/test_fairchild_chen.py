import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.ICtCp()
    # cs = colorio.cs.CAM16UCS(0.69, 20, 4.074)
    cs = colorio.cs.JzAzBz()
    # cs = colorio.cs.XYY(1)
    colorio.data.FairchildChen("SL2").show(cs)


def test_stress():
    cs = colorio.cs.CIELAB()
    ref = 4.672891455041493
    res = colorio.data.FairchildChen("SL1").stress(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref

    ref = 11.266510241793169
    res = colorio.data.FairchildChen("SL2").stress(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
    # test_stress()
