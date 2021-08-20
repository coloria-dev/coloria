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
        M_cat02 = np.array(
            [
                [+0.7328, +0.4296, -0.1624],
                [-0.7036, +1.6975, +0.0061],
                [+0.0030, +0.0136, +0.9834],
            ]
        )
        M_cat02_inv = np.linalg.inv(M_cat02)
        M_hpe = np.array(
            [
                [+0.38971, +0.68898, -0.07868],
                [-0.22981, +1.18340, +0.04641],
                [+0.00000, +0.00000, +1.00000],
            ]
        )
        M_hpe_inv = np.linalg.inv(M_hpe)
        rgb_w = M_cat02 @ self.whitepoint_xyz100

        A = np.array(
            [[0.0, 1.0, 0.0], [125 / 29, -125 / 29, 0.0], [0.0, 50 / 29, -50 / 29]]
        )
        # Ainv = np.array(
        #     [[1.0, 29 / 125, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, -116 / 200]]
        # )

        self.B = M_hpe @ M_cat02_inv @ (M_cat02.T / rgb_w).T
        self.C = A @ M_hpe_inv

        self.Binv = np.linalg.inv(self.B)
        self.Cinv = np.linalg.inv(self.C)

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        return npx.dot(self.C, f(npx.dot(self.B, xyz)))

    def to_xyz100(self, lab: ArrayLike) -> np.ndarray:
        return npx.dot(self.Binv, finv(npx.dot(self.Cinv, lab)))
