import numpy as np
from numpy.typing import ArrayLike

M16 = np.array(
    [
        [+0.401288, +0.650173, -0.051461],
        [-0.250268, +1.204414, +0.045854],
        [-0.002079, +0.048952, +0.953127],
    ]
)


def cat16(
    whitepoint_source: ArrayLike,
    whitepoint_target: ArrayLike,
    F: float,
    L_A: float,
    include_back_transform: bool = True,
    exact_inversion: bool = True,
):
    """Chromatic adaptation transform for CAM16."""
    D = F * (1 - 1 / 3.6 * np.exp((-L_A - 42) / 92))
    D = np.clip(D, 0.0, 1.0)

    whitepoint_source = np.asarray(whitepoint_source)
    whitepoint_target = np.asarray(whitepoint_target)

    rgb_w = M16 @ whitepoint_source
    Y_w = whitepoint_source[1]

    rgb_wr = M16 @ whitepoint_target
    Y_wr = whitepoint_target[1]

    D_RGB = D * (Y_w * rgb_wr) / (Y_wr * rgb_w) + 1 - D
    if exact_inversion:
        if include_back_transform:
            M = np.linalg.solve(M16, (M16.T * D_RGB).T)
            Minv = np.linalg.solve(M16, (M16.T / D_RGB).T)
        else:
            # Skip the transformation back into XYZ space. Used in CAM16.
            M = (M16.T * D_RGB).T
            Minv = np.linalg.solve(M16, np.diag(1.0 / D_RGB))
    else:
        # The CAM16 standard actually recommends using this approximation as
        # inversion operation.
        approx_inv_M16 = np.array(
            [
                [+1.86206786, -1.01125463, +0.14918677],
                [+0.38752654, +0.62144744, -0.00897398],
                [-0.01584150, -0.03412294, +1.04996444],
            ]
        )
        if include_back_transform:
            M = approx_inv_M16 @ (M16.T * D_RGB).T
            Minv = approx_inv_M16 @ (M16.T / D_RGB).T
        else:
            # Skip the transformation back into XYZ space. Used in CAM16.
            M = (M16.T * D_RGB).T
            Minv = approx_inv_M16 @ np.diag(1.0 / D_RGB)

    return M, Minv
