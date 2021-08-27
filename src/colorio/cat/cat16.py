import npx
import numpy as np
from numpy.typing import ArrayLike

M16 = np.array(
    [
        [+0.401288, +0.650173, -0.051461],
        [-0.250268, +1.204414, +0.045854],
        [-0.002079, +0.048952, +0.953127],
    ]
)


class CAT16:
    """Chromatic adaptation transform for CAM16."""

    def __init__(
        self,
        whitepoint_test: ArrayLike,
        whitepoint_reference: ArrayLike,
        F: float,
        L_A: float,
        back_transform: bool = True,
    ):
        D = F * (1 - 1 / 3.6 * np.exp((-L_A - 42) / 92))
        D = np.clip(D, 0.0, 1.0)

        rgb_w = M16 @ whitepoint_test
        Y_w = whitepoint_test[1]

        rgb_wr = M16 @ whitepoint_reference
        Y_wr = whitepoint_reference[1]

        D_RGB = D * (Y_w * rgb_wr) / (Y_wr * rgb_w) + 1 - D
        if back_transform:
            self.M = np.linalg.solve(M16, (M16.T * D_RGB).T)
        else:
            # Skip the transformation back into XYZ space. Used in CAM16.
            self.M = (M16.T * D_RGB).T

        self.Minv = np.linalg.inv(self.M)

        # The CAM16 standard actually recommends using this approximation as inversion
        # operation.
        # approx_inv_M16 = np.array(
        #     [
        #         [+1.86206786, -1.01125463, +0.14918677],
        #         [+0.38752654, +0.62144744, -0.00897398],
        #         [-0.01584150, -0.03412294, +1.04996444],
        #     ]
        # )
        # self.Minv = np.linalg.inv(self.M) if exact_inversion else approx_inv_M16 / D_RGB

    def apply(self, xyz):
        return npx.dot(self.M, xyz)

    def apply_inv(self, xyz):
        return npx.dot(self.Minv, xyz)
