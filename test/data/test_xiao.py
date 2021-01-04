import numpy
import pytest

import colorio


def test_show():
    # cs = colorio.CIELAB()
    # cs = colorio.CIEHCL()
    # cs = colorio.CIELCH()
    cs = colorio.OsaUcs()
    # cs = colorio.IPT()
    # cs = colorio.XYY()
    colorio.data.xiao.show(cs)
    # colorio.data.ebner_fairchild.savefig(cs)


@pytest.mark.parametrize(
    "cs,ref",
    [
        (colorio.CAM16UCS(0.69, 20, 64 / numpy.pi / 5), 0.7186726460882681),
        (colorio.CAM02("UCS", 0.69, 20, 64 / numpy.pi / 5), 0.7799631516138833),
        (colorio.CIEHCL(), 2.6848495918821063),
        (colorio.CIELAB(), 0.5528109307168007),
        (colorio.CIELCH(), 1.1028771830194923),
        (colorio.CIELUV(), 0.6518104837712346),
        (colorio.IPT(), 0.8928735078622902),
        (colorio.JzAzBz(), 0.7538224281799465),
        (colorio.OKLAB(), 0.758035735312263),
        (colorio.OsaUcs(), 0.8347566618822679),
        (colorio.RLAB(), 0.7089776800865354),
        (colorio.XYY(), 0.4410546230668222),
    ],
)
def test_residuals(cs, ref):
    print(cs)
    res = sum(colorio.data.xiao.residuals(cs))
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
