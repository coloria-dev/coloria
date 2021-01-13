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
    cs = colorio.cs.XYY1()
    colorio.data.luo_rigg.show(cs)
    with tempfile.TemporaryDirectory() as tmpdir:
        colorio.data.luo_rigg.savefig(Path(tmpdir) / "out.png", cs, 50)


def test_residuals():
    cs = colorio.cs.CIELAB()
    ref = 67.83352912024061
    res = colorio.data.luo_rigg.stress(cs, 8)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    test_show()
    # test_residuals()
