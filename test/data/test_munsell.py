import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    cs = colorio.cs.IPT()
    colorio.data.munsell.show(cs, V=5)
    # colorio.data.ebner_fairchild.savefig(cs)


if __name__ == "__main__":
    test_show()
