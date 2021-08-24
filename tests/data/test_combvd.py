import colorio


def test_show():
    cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.OKLAB()
    cs = colorio.cs.CAM02("UCS", 0.69, 20, 4.074)
    # cs = colorio.cs.CAM16UCS(0.69, 20, 4.074)
    # cs = colorio.cs.JzAzBz()
    # cs = colorio.cs.XYY(1)
    plt = colorio.data.COMBVD().plot(cs)
    plt.show()
    plt.close()


def test_stress():
    data = colorio.data.COMBVD()
    # cs = colorio.cs.CIELAB(data.whitepoint_xyz100)
    cs = colorio.cs.CIELAB()
    ref = 43.87405562383293
    res = data.stress(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


def test_stress_relative():
    data = colorio.data.COMBVD()
    # cs = colorio.cs.CIELAB(data.whitepoint_xyz100)
    cs = colorio.cs.CIELAB()

    ref = 44.33386329705589
    res = data.stress(cs, variant="relative")
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    # test_show()
    test_stress()
