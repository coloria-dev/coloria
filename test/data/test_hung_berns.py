import colorio


def test_show():
    # cs = colorio.CIELAB()
    # cs = colorio.CIEHCL()
    # cs = colorio.CIELCH()
    # cs = colorio.OsaUcs()
    cs = colorio.IPT()
    colorio.data.hung_berns.show(cs)
    # colorio.data.ebner_fairchild.savefig(cs)


def test_residuals():
    cs = colorio.CIELAB()
    ref = 1.178429142258768
    res = sum(colorio.data.hung_berns.residuals(cs))
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
