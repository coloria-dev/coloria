import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.OKLAB()
    cs = colorio.cs.CAM02("UCS", 0.69, 20, 4.074)
    # cs = colorio.cs.CAM16UCS(0.69, 20, 4.074)
    # cs = colorio.cs.JzAzBz()
    # cs = colorio.cs.XYY(1)
    plt = colorio.data.BfdP().plot(cs)
    plt.show()
    plt.close()


def test_stress():
    cs = colorio.cs.CIELAB()
    ref = 42.4095489785527
    res = colorio.data.BfdP().stress(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref

    ref = 44.155521908075116
    res = colorio.data.BfdP().stress(cs, variant="relative")
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    # test_show()
    test_stress()
