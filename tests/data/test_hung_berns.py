import numpy as np

import colorio


def test_show():
    cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    colorio.data.HungBerns().plot(cs).show()


def test_residuals():
    data = colorio.data.HungBerns()
    cs = colorio.cs.CIELAB(data.whitepoint_xyz100)
    ref = 4.1482146919285565
    res = np.average(data.stress(cs))
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
