import numpy as np
import pytest

import colorio


def test_show():
    cs = colorio.cs.IPT
    plt = colorio.data.EbnerFairchild().plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.OsaUcs, 3.812972375907544),
        (colorio.cs.IPT, 3.8585452834894975),
        (colorio.cs.OKLAB, 4.057308865601472),
        (colorio.cs.JzAzBz, 4.114638441215429),
        (colorio.cs.SRLAB2, 4.472312427743778),
        (colorio.cs.XYY100, 4.750458189192445),
        (colorio.cs.CIELUV, 4.813605849721629),
        (colorio.cs.CAM02LCD, 4.856397697283365),
        (colorio.cs.PROLAB, 4.932645332164939),
        (colorio.cs.CAM02UCS, 4.968202307667721),
        (colorio.cs.CAM02SCD, 5.027025387481754),
        (colorio.cs.CAM16UCS, 5.258885174710553),
        (colorio.cs.CIELAB, 5.292068327281036),
        (colorio.cs.RLAB, 5.3335671752151335),
        (colorio.cs.ICtCp, 7.048597385158547),
        (colorio.cs.CIELCH, 11.89507120054368),
        (colorio.cs.CIEHCL, 13.037122547036798),
    ],
)
def test_stress(cs_class, ref):
    print(cs_class)
    res = np.average(colorio.data.EbnerFairchild().stress(cs_class))
    print(res)
    assert abs(res - ref) < 1.0e-13 * ref


if __name__ == "__main__":
    test_show()
