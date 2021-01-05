import colorio


def test_show():
    cs = colorio.CIELAB()
    # cs = colorio.CIEHCL()
    # cs = colorio.CIELCH()
    # cs = colorio.OsaUcs()
    # cs = colorio.IPT()
    # cs = colorio.OKLAB()
    # cs = colorio.CAM02("UCS", 0.69, 20, 4.074)
    # cs = colorio.CAM16UCS(0.69, 20, 4.074)
    # cs = colorio.JzAzBz()  # TODO
    colorio.data.macadam_1942.show(cs, 50)
    colorio.data.macadam_1942.savefig("out.png", cs, 50)


def test_residuals():
    cs = colorio.CIELAB()
    ref = 0.4140811777348482
    res = colorio.data.macadam_1942.residuals(cs, 0.5)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
    # test_residuals()
