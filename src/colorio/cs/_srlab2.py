"""
SRLAB2 is an alternative to CIELAB. It simply replaces the white-point adjustment with
the more modern CIECAM02 variant. See <https://www.magnetkern.de/srlab2.html> for
details.
"""
import npx
import numpy as np
from numpy.typing import ArrayLike

from ..cat.cat02 import M_cat02
from ..illuminants import whitepoints_cie1931
from ._ciecam02 import M_hpe
from ._cielab import A, f, finv
from ._color_space import ColorSpace
from ._helpers import register


class SRLAB2(ColorSpace):
    name = "SRLAB2"
    labels = ("L", "a", "b")
    k0 = 0

    def __init__(self, whitepoint: ArrayLike = whitepoints_cie1931["D65"]):
        self.whitepoint_xyz100 = np.asarray(whitepoint)

        rgb_w = M_cat02 @ self.whitepoint_xyz100

        self.B = M_hpe @ np.linalg.inv(M_cat02) @ (M_cat02.T / rgb_w).T
        self.C = A @ np.linalg.inv(M_hpe)

        self.Binv = np.linalg.inv(self.B)
        self.Cinv = np.linalg.inv(self.C)

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        return npx.dot(self.C, f(npx.dot(self.B, xyz)))

    def to_xyz100(self, lab: ArrayLike) -> np.ndarray:
        return npx.dot(self.Binv, finv(npx.dot(self.Cinv, lab)))


register("srlab2", SRLAB2())
