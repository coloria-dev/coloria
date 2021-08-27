"""
Changjun Li, M. Ronnier Luo, Bryan Rigg, Robert W. G. Hunt,
CMC 2000 Chromatic Adaptation Transform: CMCCAT2000,
<https://doi.org/10.1002/col.10005>.
"""
import npx
import numpy as np
from numpy.typing import ArrayLike


class CMCCAT2000:
    def __init__(
        self,
        F: float,
        L_A1: float,
        L_A2: float,
        whitepoint_test: ArrayLike,
        whitepoint_reference: ArrayLike,
    ):
        """
        F=1 for average viewing condition, F=0.8 for dim- and dark-surround (projected
        image) conditions, and LA1 and LA2 are the luminances of the test and reference
        adapting fields, respectively.
        """
        M = np.array(
            [
                [0.7982, 0.3389, -0.1371],
                [-0.5918, 1.5512, 0.0406],
                [0.0008, 0.0239, 0.9753],
            ]
        )
        # Minv_approx = np.array(
        #     [
        #         [1.076450, -0.237662, 0.161212],
        #         [0.410964, 0.554342, 0.034694],
        #         [-0.010954, -0.013389, 1.024343],
        #     ]
        # )

        L_A_avg = 0.5 * (L_A1 + L_A2)
        D = F * (0.08 * np.log10(L_A_avg) + 0.76 - 0.45 * (L_A1 - L_A2) / (L_A1 + L_A2))
        D = np.clip(D, 0.0, 1.0)

        rgb_w = M @ whitepoint_test
        rgb_wr = M @ whitepoint_reference
        alpha = D * whitepoint_test[1] / whitepoint_reference[1]
        d_rgb = alpha * (rgb_wr / rgb_w) + 1 - D

        self.A = npx.solve(M, (M.T * d_rgb).T)
        self.Ainv = np.linalg.inv(self.A)

    def apply(self, xyz: ArrayLike) -> np.ndarray:
        return npx.dot(self.A, xyz)

    def apply_inv(self, xyz: ArrayLike) -> np.ndarray:
        return npx.dot(self.Ainv, xyz)
