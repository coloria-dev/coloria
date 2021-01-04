import numpy
import pytest

import colorio


def test_show():
    # cs = colorio.CIELAB()
    # cs = colorio.CIEHCL()
    # cs = colorio.CIELCH()
    # cs = colorio.OsaUcs()
    cs = colorio.IPT()
    colorio.data.hung_berns.show(cs)
    # colorio.data.ebner_fairchild.savefig(cs)


@pytest.mark.parametrize(
    "cs,ref",
    [
        (colorio.CAM16UCS(0.69, 20, 64 / numpy.pi / 5), 1.0202413285754055),
        (colorio.CAM02("UCS", 0.69, 20, 64 / numpy.pi / 5), 0.8264301197905258),
        (colorio.CIEHCL(), 3.8756472904166714),
        (colorio.CIELAB(), 1.178429142258768),
        (colorio.CIELCH(), 3.4779886316396875),
        (colorio.CIELUV(), 1.11324881108141),
        (colorio.IPT(), 0.7792185503616481),
        (colorio.JzAzBz(), 0.7079774557064785),
        (colorio.OKLAB(), 0.6989493527233173),
        (colorio.OsaUcs(), 0.6657006649233358),
        (colorio.RLAB(), 1.1409653782820464),
        (colorio.XYY(), 0.8966486255139057),
    ],
)
def test_residuals(cs, ref):
    print(cs)
    res = sum(colorio.data.hung_berns.residuals(cs))
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
