import tempfile
from pathlib import Path

import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.OKLAB()
    # cs = colorio.cs.CAM02("UCS", 0.69, 20, 4.074)
    # cs = colorio.cs.CAM16UCS(0.69, 20, 4.074)
    # cs = colorio.cs.JzAzBz()
    # cs = colorio.cs.XYY(1)
    colorio.data.MacAdam1942().show(cs)
    with tempfile.TemporaryDirectory() as tmpdir:
        colorio.data.MacAdam1942().savefig(Path(tmpdir) / "out.png", cs, 50)


def test_residuals():
    cs = colorio.cs.CIELAB()
    ref = 44.89521555157901
    res = colorio.data.MacAdam1942().stress(cs, 50.0)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
    # test_residuals()
