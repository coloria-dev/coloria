import npx
import numpy as np
from numpy.typing import ArrayLike

from ..illuminants import whitepoints_cie1931
from ._cielab import f, finv
from ._color_space import ColorSpace


class SRLAB2(ColorSpace):
    """
    SRLAB2 is an alternative to CIELAB. It simply replaces the white-point adjustment
    with the more modern CIECAM02 variant. See <https://www.magnetkern.de/srlab2.html>
    for details.
    """

    def __init__(self, whitepoint: ArrayLike = whitepoints_cie1931["D65"]):
        super().__init__("SRLAB2", ("L", "a", "b"), 0)
        self.whitepoint_xyz100 = np.asarray(whitepoint)
        self.M_cat02 = np.array(
            [
                [+0.7328, +0.4296, -0.1624],
                [-0.7036, +1.6975, +0.0061],
                [+0.0030, +0.0136, +0.9834],
            ]
        )
        self.M_cat02_inv = np.linalg.inv(self.M_cat02)
        self.M_hpe = np.array(
            [
                [+0.38971, +0.68898, -0.07868],
                [-0.22981, +1.18340, +0.04641],
                [+0.00000, +0.00000, +1.00000],
            ]
        )
        self.M_hpe_inv = np.linalg.inv(self.M_hpe)
        self.rgb_w = self.M_cat02 @ self.whitepoint_xyz100

        self.A = np.array(
            [[0.0, 1.0, 0.0], [125 / 29, -125 / 29, 0.0], [0.0, 50 / 29, -50 / 29]]
        )
        self.Ainv = np.array(
            [[1.0, 29 / 125, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, -116 / 200]]
        )

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        xyz = np.asarray(xyz)
        xyz_ = npx.dot(self.M_cat02_inv, (npx.dot(self.M_cat02, xyz).T / self.rgb_w).T)
        return npx.dot(self.A, npx.dot(self.M_hpe_inv, f(npx.dot(self.M_hpe, xyz_))))

    def to_xyz100(self, lab: ArrayLike) -> np.ndarray:
        lab = np.asarray(lab)
        xyz_ = npx.dot(
            self.M_hpe_inv, finv(npx.dot(self.M_hpe, npx.dot(self.Ainv, lab)))
        )
        return npx.dot(self.M_cat02_inv, (npx.dot(self.M_cat02, xyz_).T * self.rgb_w).T)
