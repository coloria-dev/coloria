import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    cs = colorio.cs.IPT()
    colorio.data.ebner_fairchild.show(cs)
    # colorio.data.ebner_fairchild.savefig(cs)


def test_residuals():
    cs = colorio.cs.CIELAB()
    ref = 4.095516235268518
    res = sum(colorio.data.ebner_fairchild.residuals(cs))
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
