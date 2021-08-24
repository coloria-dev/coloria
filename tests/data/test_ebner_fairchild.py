import numpy as np

import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    cs = colorio.cs.IPT()
    plt = colorio.data.EbnerFairchild().plot(cs)
    plt.show()
    plt.close()


def test_stress():
    cs = colorio.cs.CIELAB()
    ref = 5.3071509533648085
    res = np.average(colorio.data.EbnerFairchild().stress(cs))
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
