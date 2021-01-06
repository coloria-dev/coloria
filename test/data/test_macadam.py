import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.OKLAB()
    # cs = colorio.cs.CAM02("UCS", 0.69, 20, 4.074)
    # cs = colorio.cs.CAM16UCS(0.69, 20, 4.074)
    # cs = colorio.cs.JzAzBz()
    cs = colorio.cs.XYY()
    colorio.data.macadam_1942.show(cs, 50)
    colorio.data.macadam_1942.savefig("out.png", cs, 50)


def test_residuals():
    cs = colorio.cs.CIELAB()
    ref = 10.418677911638703
    res = colorio.data.macadam_1942.residuals(cs, 0.5)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
    # test_residuals()
