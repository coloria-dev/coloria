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
    def __init__(self, whitepoint_test: ArrayLike, whitepoint_reference: ArrayLike):
        judd = np.array(
            [
                [0.0, 1.0, 0.0],
                [-0.46, 1.36, 0.1],
                [0.0, 0.0, 1.0],
            ]
        )
        # judd_inv_approx = np.array(
        #     [
        #         [2.954, -2.174, 0.22],
        #         [1.0, 0.0, 0.0],
        #         [0.0, 0.0, 1.0],
        #     ]
        # )
        # judd_inv = np.linalg.inv(judd) if exact_inversion else judd_inv_approx

        rgb_rw = judd @ whitepoint_reference
        rgb_w = judd @ whitepoint_test
        abc = rgb_rw / rgb_w

        self.M = np.linalg.solve(judd, (judd.T * abc).T)
        self.Minv = np.linalg.inv(self.M)

    def apply(self, xyz: ArrayLike) -> np.ndarray:
        return npx.dot(self.M, xyz)

    def apply_inv(self, xyz: ArrayLike) -> np.ndarray:
        return npx.dot(self.Minv, xyz)
