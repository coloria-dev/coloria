"""
See

Ming Ronnier Luo,
CIE Chromatic Adaptation; Comparison of von Kries, CIELAB, CMCCAT97 and CAT02,
<https://doi.org/10.1007/978-3-642-27851-8_321-1>.
"""
import numpy as np
from numpy.typing import ArrayLike


class CMCCAT97:
    def __init__(
        self,
        F: float,
        L_A: float,
        whitepoint_test: ArrayLike,
        whitepoint_reference: ArrayLike,
    ):
        self.M = np.array(
            [
                [0.8951, 0.2664, -0.1614],
                [-0.7502, 1.7135, 0.0367],
                [0.0389, -0.0685, 1.0296],
            ]
        )
        self.Minv = np.linalg.inv(self.M)

        D = F - F / (1 + 2 * L_A ** 0.25 + L_A ** 2 / 300)

        r_w, g_w, b_w = self.M @ whitepoint_test
        r_rw, g_rw, b_rw = self.M @ whitepoint_reference

        self.p = (b_w / b_rw) ** 0.0834
        ratio = np.array([r_rw / r_w, g_rw / g_w, b_rw / b_w ** self.p])
        self.d_rgb = D * ratio + 1 - D

    def apply(self, xyz: ArrayLike) -> np.ndarray:
        xyz = np.asarray(xyz)
        Y = xyz[1].copy()
        rgb = self.M @ (xyz / Y)
        rgb[2] = np.sign(rgb[2]) * np.abs(rgb[2]) ** self.p
        rgb_c = self.d_rgb * rgb
        return self.Minv @ (rgb_c * Y)

    # def apply_inv(self, xyz: ArrayLike) -> np.ndarray:
    #     return npx.dot(self.judd_inv, (np.dot(self.judd, xyz).T / self.abc).T)
