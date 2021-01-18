import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    cs = colorio.cs.IPT()
    colorio.data.Munsell().show(cs, V=5)


def test_stress():
    cs = colorio.cs.CIELAB()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.JzAzBz()
    ref = 0.6928265531073651
    res = colorio.data.Munsell().stress_lightness(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_stress()
