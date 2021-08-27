"""
See

Ming Ronnier Luo,
CIE Chromatic Adaptation; Comparison of von Kries, CIELAB, CMCCAT97 and CAT02,
<https://doi.org/10.1007/978-3-642-27851-8_321-1>.
"""
import npx
import numpy as np
from numpy.typing import ArrayLike


class VonKries:
    def __init__(
        self,
        whitepoint_test: ArrayLike,
        whitepoint_reference: ArrayLike,
        exact_inversion: bool = True,
    ):
        self.judd = np.array([[0.0, 1.0, 0.0], [-0.46, 1.36, 0.1], [0.0, 0.0, 1.0]])
        judd_inv_approx = np.array(
            [[2.954, -2.174, 0.22], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        )
        self.judd_inv = np.linalg.inv(self.judd) if exact_inversion else judd_inv_approx

        rgb_rw = self.judd @ whitepoint_reference
        rgb_w = self.judd @ whitepoint_test
        self.abc = rgb_rw / rgb_w

    def apply(self, xyz: ArrayLike) -> np.ndarray:
        return npx.dot(self.judd_inv, (np.dot(self.judd, xyz).T * self.abc).T)

    def apply_inv(self, xyz: ArrayLike) -> np.ndarray:
        return npx.dot(self.judd_inv, (np.dot(self.judd, xyz).T / self.abc).T)
