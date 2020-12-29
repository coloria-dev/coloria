import colorio



def test_ebner_fairchild():
    cs = colorio.CIELAB()
    colorio.data.ebner_fairchild.show(cs)
    # colorio.data.ebner_fairchild.savefig(cs)
    # colorio.data.ebner_fairchild.residual(cs)


if __name__ == "__main__":
    test_ebner_fairchild()
