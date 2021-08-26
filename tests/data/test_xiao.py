import numpy as np
import pytest

import colorio


def test_show():
    cs = colorio.cs.OsaUcs
    plt = colorio.data.Xiao().plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref1,ref2",
    [
        (colorio.cs.XYY100, 3.407521864323334, 5.263310866198113),
        (colorio.cs.ICtCp, 4.167052243100838, 5.355499404429293),
        (colorio.cs.RLAB, 4.87775326766349, 6.641721031205064),
        (colorio.cs.CIELAB, 5.307352387184075, 6.544002194941983),
        (colorio.cs.CIELUV, 5.032609665293969, 7.496109721292327),
        (colorio.cs.PROLAB, 5.10727470695892, 7.486202482576113),
        (colorio.cs.JzAzBz, 5.793500094583313, 9.291838929886403),
        (colorio.cs.OKLAB, 5.832063085506336, 7.467968194998472),
        (colorio.cs.SRLAB2, 6.114067236635034, 8.417170601378004),
        (colorio.cs.CAM02LCD, 6.275569700740668, 12.226442745598048),
        (colorio.cs.OsaUcs, 6.378661932445455, 10.314700143107958),
        (colorio.cs.CAM02UCS, 6.581926127571111, 12.720448217181325),
        (colorio.cs.CAM02SCD, 6.730819889889601, 12.932583590655142),
        (colorio.cs.CAM16UCS, 6.569588788212987, 13.140364725832729),
        (colorio.cs.IPT, 6.77347388639466, 13.387243439462512),
        # (colorio.cs.CIELCH, 18.48385673790004, 31.190988974730637),
        (colorio.cs.CIEHCL, 19.709578352631745, 25.628130654353743),
    ],
)
def test_stress(cs_class, ref1, ref2):
    print(cs_class)
    vals = colorio.data.Xiao().stress(cs_class)
    res = np.average(vals)
    print(res)
    assert abs(res - ref1) < 1.0e-12 * ref1

    res = np.max(vals)
    print(res)
    assert abs(res - ref2) < 1.0e-12 * ref2


if __name__ == "__main__":
    test_show()
