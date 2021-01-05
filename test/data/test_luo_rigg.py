import colorio


def test_show():
    cs = colorio.CIELAB()
    # cs = colorio.CIEHCL()
    # cs = colorio.CIELCH()
    # cs = colorio.OsaUcs()
    # cs = colorio.IPT()
    # cs = colorio.CAM16UCS(0.69, 20, 4.074)
    colorio.data.luo_rigg.show(cs, 50)
    colorio.data.luo_rigg.savefig("out.png", cs, 50)


def test_residuals():
    cs = colorio.CIELAB()
    ref = 0.4140811777348482
    res = colorio.data.luo_rigg.residuals(cs, 0.5)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
    # test_residuals()
