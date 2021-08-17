import npx
import numpy as np
from numpy.typing import ArrayLike

from ..illuminants import whitepoints_cie1931
from ._color_space import ColorSpace


def f(t):
    delta = 6 / 29
    out = np.array(t, dtype=float)
    is_greater = out > delta ** 3
    out[is_greater] = 116 * np.cbrt(out[is_greater]) - 16
    out[~is_greater] = out[~is_greater] / (delta / 2) ** 3
    return out


def finv(t):
    t /= 116
    delta = 6 / 29
    out = np.array(t, dtype=float)
    is_greater = out > 2 / 29
    out[is_greater] = (out[is_greater] + 4 / 29) ** 3
    out[~is_greater] = 3 * delta ** 2 * out[~is_greater]
    return out


class CIELAB(ColorSpace):
    def __init__(self, whitepoint: ArrayLike = whitepoints_cie1931["D65"]):
        super().__init__("CIELAB", ("L*", "a*", "b*"), 0)
        self.whitepoint_xyz100 = np.asarray(whitepoint)
        self.whitepoint = np.array([100.0, 0.0, 0.0])

        # See
        # https://gist.github.com/nschloe/ad1d288917a140978f7db6c401cb7f17
        # for a speed comparison. None is really faster, but the matrix-approach is
        # easiest to read.
        self.M = np.array(
            [[0.0, 1.0, 0.0], [125 / 29, -125 / 29, 0.0], [0.0, 50 / 29, -50 / 29]]
        )
        self.Minv = np.array(
            [[1.0, 29 / 125, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, -116 / 200]]
        )

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        xyz = np.asarray(xyz)
        return npx.dot(self.M, f((xyz.T / self.whitepoint_xyz100).T))

    def to_xyz100(self, lab: ArrayLike) -> np.ndarray:
        lab = np.asarray(lab)
        return (finv(npx.dot(self.Minv, lab)).T * self.whitepoint_xyz100).T
