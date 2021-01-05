import colorio


def test_show():
    cs = colorio.CIELAB()
    # cs = colorio.CIEHCL()
    # cs = colorio.CIELCH()
    # cs = colorio.OsaUcs()
    # cs = colorio.IPT()
    colorio.data.macadam_1942.show(cs, 50)
    # colorio.data.macadam.savefig(cs)


def test_residuals():
    cs = colorio.CIELAB()
    ref = 0.4140811777348482
    res = colorio.data.macadam_1942.residuals(cs, 0.5)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    # test_show()
    test_residuals()
