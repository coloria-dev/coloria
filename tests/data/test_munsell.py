import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    cs = colorio.cs.IPT()
    plt = colorio.data.Munsell().plot(cs, V=5)
    plt.show()
    plt.close()


def test_show_lightness():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    cs = colorio.cs.IPT()
    # cs = colorio.cs.OKLAB()
    plt = colorio.data.Munsell().plot_lightness(cs)
    plt.show()
    plt.close()


def test_stress():
    cs = colorio.cs.CIELAB()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.JzAzBz()
    ref = 0.6928265531073651
    res = colorio.data.Munsell().stress_lightness(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    # test_stress()
    # test_show()
    test_show_lightness()
