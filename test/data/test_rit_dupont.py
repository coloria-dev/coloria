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
    # cs = colorio.cs.XYY(1)
    colorio.data.rit_dupont.show(cs, "light brown")
    with tempfile.TemporaryDirectory() as tmpdir:
        colorio.data.rit_dupont.savefig(Path(tmpdir) / "out.png", cs, "blue")


def test_residual():
    cs = colorio.cs.CIELAB()
    ref = 32.85162431714214
    res = colorio.data.rit_dupont.stress(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
    # test_residual()
