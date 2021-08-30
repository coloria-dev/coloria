import npx
import numpy as np
from numpy.typing import ArrayLike

M_cat02 = np.array(
    [
        [+0.7328, +0.4296, -0.1624],
        [-0.7036, +1.6975, +0.0061],
        [+0.0030, +0.0136, +0.9834],
    ]
)


class CAT02:
    """Chromatic adaptation transform for CIECAM02."""

    def __init__(
        self,
        whitepoint_test: ArrayLike,
        whitepoint_reference: ArrayLike,
        F: float,
        L_A: float,
    ):
        D = F * (1.0 - np.exp((-L_A - 42) / 92) / 3.6)
        D = np.clip(D, 0.0, 1.0)

        rgb_wr = M_cat02 @ whitepoint_reference
        rgb_w = M_cat02 @ whitepoint_test
        Y_w = whitepoint_test[1]
        Y_wr = whitepoint_reference[1]

        D_RGB = D * (Y_w * rgb_wr) / (Y_wr * rgb_w) + 1 - D

        M_cat02_inv = np.linalg.inv(M_cat02)
        self.M = M_cat02_inv @ (M_cat02.T * D_RGB).T
        self.Minv = M_cat02_inv @ (M_cat02.T / D_RGB).T

    def apply(self, xyz):
        return npx.dot(self.M, xyz)

    def apply_inv(self, xyz):
        return npx.dot(self.Minv, xyz)
