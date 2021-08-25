import numpy as np
import pytest

import colorio


def test_show():
    cs = colorio.cs.CIELAB
    plt = colorio.data.HungBerns().plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.OsaUcs, 2.4953099801314433),
        (colorio.cs.OKLAB, 2.633719505192444),
        (colorio.cs.JzAzBz, 2.685364268944074),
        (colorio.cs.IPT, 2.9211968310535603),
        (colorio.cs.SRLAB2, 2.9922958127771158),
        (colorio.cs.CAM02LCD, 3.047547351029614),
        (colorio.cs.CAM02UCS, 3.188456200080372),
        (colorio.cs.CAM02SCD, 3.2570190989569334),
        (colorio.cs.XYY100, 3.4401874865344246),
        (colorio.cs.CAM16UCS, 3.5285507819152624),
        (colorio.cs.CIELUV, 3.69928776179375),
        (colorio.cs.PROLAB, 3.7398740774298758),
        (colorio.cs.ICtCp, 3.7957739893557014),
        (colorio.cs.RLAB, 3.883688326905166),
        (colorio.cs.CIELAB, 4.1482146919285565),
        (colorio.cs.CIELCH, 13.833274990951756),
        (colorio.cs.CIEHCL, 14.703269694844176),
    ],
)
def test_stress(cs_class, ref):
    print(cs_class)
    res = np.average(colorio.data.HungBerns().stress(cs_class))
    print(res)
    assert abs(res - ref) < 1.0e-13 * ref


if __name__ == "__main__":
    test_show()
