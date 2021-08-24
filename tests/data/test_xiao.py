import numpy as np
import pytest

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


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.CIELUV, 2.381285478701682),
        (colorio.cs.XYY100, 2.3985002202066097),
        (colorio.cs.CIELAB, 2.788362604236149),
        (colorio.cs.PROLAB, 2.8038926193818785),
        (colorio.cs.SRLAB2, 3.1075993110078013),
        (colorio.cs.RLAB, 3.6814827262133005),
        (colorio.cs.ICtCp, 3.924094855240532),
        (colorio.cs.CAM02SCD, 4.556735406166622),
        (colorio.cs.CAM02UCS, 4.58825706817219),
        (colorio.cs.CAM02LCD, 4.618341690395547),
        (colorio.cs.CAM16UCS, 4.626957018528581),
        (colorio.cs.JzAzBz, 5.230529469594494),
        (colorio.cs.OKLAB, 6.256386444491906),
        (colorio.cs.OsaUcs, 8.22761086694146),
        (colorio.cs.IPT, 8.853878081896712),
        (colorio.cs.CIELCH, 9.983951212967135),
        (colorio.cs.CIEHCL, 11.516088298979938),
    ],
)
def test_stress(cs_class, ref):
    res = np.average(colorio.data.Xiao().stress(cs_class))
    print(cs_class)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
