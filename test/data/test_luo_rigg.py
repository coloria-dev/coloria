import tempfile
from pathlib import Path

import colorio


def test_show():
    # cs = colorio.cs.CIELAB()
    # cs = colorio.cs.CIEHCL()
    # cs = colorio.cs.CIELCH()
    # cs = colorio.cs.OsaUcs()
    # cs = colorio.cs.IPT()
    # cs = colorio.cs.CAM16UCS(0.69, 20, 4.074)
    cs = colorio.cs.XYY(1)
    colorio.data.LuoRigg(8).show(cs)
    with tempfile.TemporaryDirectory() as tmpdir:
        colorio.data.LuoRigg(8).savefig(Path(tmpdir) / "out.png", cs)


def test_residuals():
    cs = colorio.cs.CIELAB()
    ref = 48.472622852849376
    res = colorio.data.LuoRigg(8).stress(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
    # test_residuals()
