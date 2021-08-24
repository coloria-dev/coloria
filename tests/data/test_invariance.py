import numpy as np
import pytest

import colorio

rng = np.random.default_rng(0)


def dot(a, b):
    """Take arrays `a` and `b` and form the dot product between the last axis
    of `a` and the first of `b`.
    """
    b = np.asarray(b)
    return np.dot(a, b.reshape(b.shape[0], -1)).reshape(a.shape[:-1] + b.shape[1:])


class CielabScaled:
    def __init__(self, whitepoint):
        self.cielab = colorio.cs.CIELAB(whitepoint)
        self.alpha = 2.41
        self.k0 = self.cielab.k0

    def from_xyz100(self, xyz100):
        return self.cielab.from_xyz100(xyz100) * self.alpha

    def to_xyz100(self, lab):
        return self.cielab.to_xyz100(lab) / self.alpha


class CielabTranslated:
    def __init__(self, whitepoint):
        self.cielab = colorio.cs.CIELAB(whitepoint)
        self.x = rng.random(3)
        self.k0 = self.cielab.k0

    def from_xyz100(self, xyz100):
        return (self.cielab.from_xyz100(xyz100).T + self.x).T

    def to_xyz100(self, xyz100):
        return self.cielab.to_xyz100(xyz100 - self.x)


class CielabRotated:
    def __init__(self, whitepoint):
        self.cielab = colorio.cs.CIELAB(whitepoint)
        self.R, _ = np.linalg.qr(rng.random((3, 3)))
        self.Rinv = np.linalg.inv(self.R)
        self.k0 = self.cielab.k0

    def from_xyz100(self, xyz100):
        return dot(self.R, self.cielab.from_xyz100(xyz100))

    def to_xyz100(self, lab):
        return self.cielab.to_xyz100(dot(self.Rinv, lab))


@pytest.mark.parametrize(
    "fun",
    [
        colorio.data.MacAdam1974().stress,
        colorio.data.Witt().stress,
    ],
)
@pytest.mark.parametrize("ct", [CielabScaled, CielabTranslated, CielabRotated])
def test_invariance(fun, ct):
    cs = colorio.cs.CIELAB

    val0 = fun(cs)
    val1 = fun(ct)
    print(fun, ct)
    print(val0)
    print(val1)
    assert np.all(np.abs(val0 - val1) < 1.0e-14 * np.abs(val0))


# Lightness and hue linearity are are not rotation invariant because there is a
# preferred direction, the lightness direction.
@pytest.mark.parametrize(
    "fun",
    [
        colorio.data.EbnerFairchild().stress,
        colorio.data.HungBerns().stress,
        colorio.data.Xiao().stress,
        #
        colorio.data.FairchildChen("SL1").stress,
        colorio.data.FairchildChen("SL2").stress,
        colorio.data.Munsell().stress_lightness,
    ],
)
@pytest.mark.parametrize("ct", [CielabScaled, CielabTranslated])
def test_lightness(fun, ct):
    cs = colorio.cs.CIELAB

    val0 = fun(cs)
    val1 = fun(ct)
    print(fun, ct)
    print(val0)
    print(val1)
    assert np.all(np.abs(val0 - val1) < 1.0e-14 * np.abs(val0))


if __name__ == "__main__":
    test_invariance(colorio.data.ebner_fairchild)
