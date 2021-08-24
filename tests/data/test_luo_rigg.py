import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.CAM16UCS(0.69, 20, 4.074)
    cs = colorio.cs.XYY(1)
    plt = colorio.data.LuoRigg(8).plot(cs)
    plt.show()
    plt.close()


def test_residuals():
    data = colorio.data.LuoRigg(8)
    # cs = colorio.cs.CIELAB(data.whitepoint_xyz100)
    cs = colorio.cs.CIELAB()
    ref = 48.472622852849376
    res = data.stress(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
    # test_residuals()
