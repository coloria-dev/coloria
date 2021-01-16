import numpy as np
import pytest

import colorio

np.random.seed(0)


class CielabScaled:
    def __init__(self):
        self.cielab = colorio.cs.CIELAB()
        self.alpha = 2.41
        self.k0 = self.cielab.k0

    def from_xyz100(self, xyz100):
        return self.cielab.from_xyz100(xyz100) * self.alpha

    def to_xyz100(self, lab):
        return self.cielab.to_xyz100(lab) / self.alpha


class CielabTranslated:
    def __init__(self):
        self.cielab = colorio.cs.CIELAB()
        self.x = np.random.rand(3)
        self.k0 = self.cielab.k0

    def from_xyz100(self, xyz100):
        return (self.cielab.from_xyz100(xyz100).T + self.x).T

    def to_xyz100(self, xyz100):
        return self.cielab.to_xyz100(xyz100 - self.x)


class CielabRotated:
    def __init__(self):
        self.cielab = colorio.cs.CIELAB()
        self.R, _ = np.linalg.qr(np.random.rand(3, 3))
        self.Rinv = np.linalg.inv(self.R)
        self.k0 = self.cielab.k0

    def from_xyz100(self, xyz100):
        return self.R @ self.cielab.from_xyz100(xyz100)

    def to_xyz100(self, lab):
        return self.cielab.to_xyz100(self.Rinv @ lab)


@pytest.mark.parametrize(
    "fun",
    [
        colorio.data.ebner_fairchild.stress,
        colorio.data.hung_berns.stress,
        colorio.data.macadam_1974.stress,
        colorio.data.witt.stress,
        colorio.data.xiao.stress,
    ],
)
@pytest.mark.parametrize("ct", [CielabScaled(), CielabTranslated(), CielabRotated()])
def test_invariance(fun, ct):
    cs = colorio.cs.CIELAB()

    val0 = fun(cs)
    val1 = fun(ct)
    print(fun, ct)
    print(val0)
    print(val1)
    assert np.all(np.abs(val0 - val1) < 1.0e-14 * np.abs(val0))


# test lightness separately because they are not rotation invariant
@pytest.mark.parametrize(
    "fun",
    [
        lambda cs: colorio.data.fairchild_chen.stress(cs, "SL1"),
        lambda cs: colorio.data.fairchild_chen.stress(cs, "SL2"),
        colorio.data.munsell.stress_lightness,
    ],
)
@pytest.mark.parametrize("ct", [CielabScaled(), CielabTranslated()])
def test_lightness(fun, ct):
    cs = colorio.cs.CIELAB()

    val0 = fun(cs)
    val1 = fun(ct)
    print(fun, ct)
    print(val0)
    print(val1)
    assert np.all(np.abs(val0 - val1) < 1.0e-14 * np.abs(val0))


if __name__ == "__main__":
    test_invariance(colorio.data.ebner_fairchild)
