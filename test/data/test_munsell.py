import colorio


def test_show():
    # cs = colorio.CIELAB()
    # cs = colorio.CIEHCL()
    # cs = colorio.CIELCH()
    # cs = colorio.OsaUcs()
    cs = colorio.IPT()
    colorio.data.munsell.show(cs, V=5)
    # colorio.data.ebner_fairchild.savefig(cs)


if __name__ == "__main__":
    test_show()
