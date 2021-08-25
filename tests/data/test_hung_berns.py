import numpy as np
import pytest

import colorio


def test_show():
    cs = colorio.cs.CIELAB
    plt = colorio.data.HungBerns().plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref1,ref2",
    [
        (colorio.cs.OsaUcs, 2.4953099801314433, 6.59231944961339),
        (colorio.cs.OKLAB, 2.633719505192444, 9.319977779587065),
        (colorio.cs.JzAzBz, 2.685364268944074, 5.55970441646704),
        (colorio.cs.IPT, 2.9211968310535603, 5.89478443004315),
        (colorio.cs.SRLAB2, 2.9922958127771158, 7.93281663197524),
        (colorio.cs.CAM02LCD, 3.047547351029614, 7.487530741518032),
        (colorio.cs.CAM02UCS, 3.188456200080372, 7.778039133853468),
        (colorio.cs.CAM02SCD, 3.2570190989569334, 7.928007054988701),
        (colorio.cs.XYY100, 3.4401874865344246, 6.187170189969933),
        (colorio.cs.CAM16UCS, 3.5285507819152624, 16.03234710518988),
        (colorio.cs.CIELUV, 3.69928776179375, 8.563919408909138),
        (colorio.cs.PROLAB, 3.7398740774298758, 9.334432349567447),
        (colorio.cs.ICtCp, 3.7957739893557014, 11.189083034653905),
        (colorio.cs.RLAB, 3.883688326905166, 13.601209583932523),
        (colorio.cs.CIELAB, 4.1482146919285565, 15.031157563615553),
        # (colorio.cs.CIELCH, 13.833274990951756, 1.0),
        (colorio.cs.CIEHCL, 14.703269694844176, 23.519088290340893),
    ],
)
def test_stress(cs_class, ref1, ref2):
    print(cs_class)
    vals = colorio.data.HungBerns().stress(cs_class)
    res = np.average(vals)
    print(res)
    assert abs(res - ref1) < 1.0e-13 * ref1

    res = np.max(vals)
    print(res)
    assert abs(res - ref2) < 1.0e-13 * ref2


if __name__ == "__main__":
    test_show()
