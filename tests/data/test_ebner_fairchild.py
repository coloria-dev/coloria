import numpy as np
import pytest

import colorio


def test_show():
    cs = colorio.cs.IPT
    plt = colorio.data.EbnerFairchild().plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref1,ref2",
    [
        (colorio.cs.OsaUcs, 3.812972375907544, 7.178611330824679),
        (colorio.cs.IPT, 3.8585452834894975, 7.722339571500826),
        (colorio.cs.OKLAB, 4.057308865601472, 7.161717923241752),
        (colorio.cs.JzAzBz, 4.114638441215429, 8.80555808428739),
        (colorio.cs.SRLAB2, 4.472312427743778, 8.872702833107542),
        (colorio.cs.XYY100, 4.750458189192445, 13.435491691097537),
        (colorio.cs.CIELUV, 4.813605849721629, 8.914484060462362),
        (colorio.cs.CAM02LCD, 4.856397697283365, 8.050603730839033),
        (colorio.cs.PROLAB, 4.932645332164939, 8.4138784020548),
        (colorio.cs.CAM02UCS, 4.968202307667721, 8.220723083292194),
        (colorio.cs.CAM02SCD, 5.027025387481754, 8.279318059320925),
        (colorio.cs.CAM16UCS, 5.258885174710553, 8.660220687471277),
        (colorio.cs.CIELAB, 5.292068327281036, 10.456987379313865),
        (colorio.cs.RLAB, 5.3335671752151335, 9.965820911777676),
        (colorio.cs.ICtCp, 7.048597385158547, 13.628048776065693),
        # (colorio.cs.CIELCH, 11.89507120054368, 1.0),
        (colorio.cs.CIEHCL, 13.037122547036798, 20.663099839073062),
    ],
)
def test_stress(cs_class, ref1, ref2):
    print(cs_class)
    vals = colorio.data.EbnerFairchild().stress(cs_class)
    res = np.average(vals)
    print(res)
    assert abs(res - ref1) < 1.0e-11 * ref1

    res = np.max(vals)
    print(res)
    assert abs(res - ref2) < 1.0e-11 * ref2


if __name__ == "__main__":
    test_show()
