import numpy

import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.JzAzBz()
    cs = colorio.cs.XYY(1)
    colorio.data.fairchild_chen.show(cs, "SL2")
    # colorio.data.fairchild_chen.savefig(cs)


def test_stress():
    cs = colorio.cs.CIELAB()
    ref = 4.672891455041493
    res = colorio.data.fairchild_chen.stress(cs, "SL1")
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref

    ref = 11.266510241793169
    res = colorio.data.fairchild_chen.stress(cs, "SL2")
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    # test_show()
    test_stress()
