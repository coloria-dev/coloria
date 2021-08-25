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
        (colorio.cs.CIELUV, 2.381285478701682, 3.1332518106519083),
        (colorio.cs.XYY100, 2.3985002202066097, 3.4704952491749737),
        (colorio.cs.CIELAB, 2.788362604236149, 5.414285724429276),
        (colorio.cs.PROLAB, 2.8038926193818785, 4.372926666259824),
        (colorio.cs.SRLAB2, 3.1075993110078013, 6.141156849717835),
        (colorio.cs.RLAB, 3.6814827262133005, 5.822479304956172),
        (colorio.cs.ICtCp, 3.924094855240532, 5.526529850337941),
        (colorio.cs.CAM02SCD, 4.556735406166622, 7.79989808564784),
        (colorio.cs.CAM02UCS, 4.58825706817219, 7.836629887470947),
        (colorio.cs.CAM02LCD, 4.618341690395547, 7.839432495812633),
        (colorio.cs.CAM16UCS, 4.626957018528581, 7.639901148051151),
        (colorio.cs.JzAzBz, 5.230529469594494, 10.855292125203963),
        (colorio.cs.OKLAB, 6.256386444491906, 9.59481014620974),
        (colorio.cs.OsaUcs, 8.22761086694146, 12.136891378166084),
        (colorio.cs.IPT, 8.853878081896712, 16.22085669819517),
        # (colorio.cs.CIELCH, 9.983951212967135, 1.0),
        (colorio.cs.CIEHCL, 11.516088298979938, 21.378029058252057),
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
