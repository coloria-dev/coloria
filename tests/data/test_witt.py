import pytest

import colorio


def test_show():
    cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.OKLAB()
    cs = colorio.cs.CAM02("UCS", 0.69, 20, 4.074)
    # cs = colorio.cs.CAM16UCS(0.69, 20, 4.074)
    # cs = colorio.cs.JzAzBz()
    # cs = colorio.cs.XYY(1)
    plt = colorio.data.Witt().plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.CIELAB, 51.70889935711488),
        (colorio.cs.CAM02UCS, 30.4584731112042),
        (colorio.cs.CAM02LCD, 39.16831128020275),
        (colorio.cs.CAM02SCD, 28.577269461047187),
        (colorio.cs.OKLAB, 45.2276225067436),
        (colorio.cs.JzAzBz, 48.53002773753702),
        (colorio.cs.CAM16UCS, 30.991751511774744),
        (colorio.cs.SRLAB2, 51.549289192015046),
    ],
)
def test_stress(cs_class, ref):
    res = colorio.data.Witt().stress(cs_class)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref
