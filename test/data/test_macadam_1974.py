import tempfile
from pathlib import Path

import colorio


def test_show():
    cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.OKLAB()
    # cs = colorio.cs.CAM02("UCS", 0.69, 20, 4.074)
    # cs = colorio.cs.CAM16UCS(0.69, 20, 4.074)
    # cs = colorio.cs.JzAzBz()
    # cs = colorio.cs.XYY1()
    colorio.data.macadam_1974.show(cs)
    with tempfile.TemporaryDirectory() as tmpdir:
        colorio.data.macadam_1974.savefig(Path(tmpdir) / "out.png", cs, 50)


def test_residual():
    cs = colorio.cs.CIELAB()
    ref = 0.2453191916738762
    res = colorio.data.macadam_1974.residuals(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
    # test_residual()
