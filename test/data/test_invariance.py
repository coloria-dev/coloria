import numpy as np
import pytest

import colorio


class CielabTranslate():
    def __init__(self):
        self.cielab = colorio.cs.CIELAB()
        self.x = np.array([0.1, 0.4, -1.1])
        self.k0 = self.cielab.k0

    def from_xyz100(self, xyz100):
        return (self.cielab.from_xyz100(xyz100).T + self.x).T

    def to_xyz100(self, xyz100):
        return self.cielab.to_xyz100(xyz100 - self.x)


@pytest.mark.parametrize("fun", [
    colorio.data.ebner_fairchild.stress,
    lambda cs: colorio.data.fairchild_chen.stress(cs, "SL1"),
    lambda cs: colorio.data.fairchild_chen.stress(cs, "SL2"),
    colorio.data.hung_berns.stress,
    colorio.data.macadam_1974.stress,
    colorio.data.witt.stress,
    colorio.data.xiao.stress,
])
def test_invariance(fun):
    cs = colorio.cs.CIELAB()
    ct = CielabTranslate()

    val0 = fun(cs)
    val1 = fun(ct)
    print(val0)
    print(val1)
    assert np.all(np.abs(val0 - val1) < 1.0e-14 * np.abs(val0))



if __name__ == "__main__":
    test_invariance(colorio.data.ebner_fairchild)
