import numpy
import pytest

import colorio


def test_show():
    # cs = colorio.CIELAB()
    # cs = colorio.CIEHCL()
    # cs = colorio.CIELCH()
    # cs = colorio.OsaUcs()
    cs = colorio.IPT()
    colorio.data.ebner_fairchild.show(cs)
    # colorio.data.ebner_fairchild.savefig(cs)


@pytest.mark.parametrize(
    "cs,ref",
    [
        (colorio.CAM16UCS(0.69, 20, 64 / numpy.pi / 5), 3.7383772675948763),
        (colorio.CAM02("UCS", 0.69, 20, 64 / numpy.pi / 5), 3.492838139468934),
        (colorio.CIEHCL(), 12.619950376639341),
        (colorio.CIELAB(), 4.095516235268518),
        (colorio.CIELCH(), 13.405076662924223),
        (colorio.CIELUV(), 3.7225918914653073),
        (colorio.IPT(), 2.983536242312675),
        (colorio.JzAzBz(), 3.1235467308302223),
        (colorio.OKLAB(), 3.0983054621799875),
        (colorio.OsaUcs(), 2.948843253469774),
        (colorio.RLAB(), 4.1766551329386346),
        (colorio.XYY(), 3.5552420588236453),
    ],
)
def test_residuals(cs, ref):
    print(cs)
    res = sum(colorio.data.ebner_fairchild.residuals(cs))
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
