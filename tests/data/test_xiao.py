import numpy as np

import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.XYY()
    plt = colorio.data.Xiao().plot(cs)
    plt.show()
    plt.close()


def test_residuals():
    data = colorio.data.Xiao()
    cs = colorio.cs.CIELAB(whitepoint=data.whitepoint_xyz100)
    ref = 2.7883626042361533
    res = np.average(data.stress(cs))
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
