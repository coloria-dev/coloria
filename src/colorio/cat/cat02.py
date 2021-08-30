import numpy as np
from numpy.typing import ArrayLike

M_cat02 = np.array(
    [
        [+0.7328, +0.4296, -0.1624],
        [-0.7036, +1.6975, +0.0061],
        [+0.0030, +0.0136, +0.9834],
    ]
)


def cat02(
    whitepoint_source: ArrayLike,
    whitepoint_target: ArrayLike,
    F: float,
    L_A: float,
    exact_inversion: bool = True,
):
    """Chromatic adaptation transform for CIECAM02."""
    D = F * (1.0 - np.exp((-L_A - 42) / 92) / 3.6)
    D = np.clip(D, 0.0, 1.0)

    whitepoint_source = np.asarray(whitepoint_source)
    whitepoint_target = np.asarray(whitepoint_target)

    rgb_wr = M_cat02 @ whitepoint_target
    rgb_w = M_cat02 @ whitepoint_source
    Y_w = whitepoint_source[1]
    Y_wr = whitepoint_target[1]

    D_RGB = D * (Y_w * rgb_wr) / (Y_wr * rgb_w) + 1 - D

    if exact_inversion:
        M = np.linalg.solve(M_cat02, (M_cat02.T * D_RGB).T)
        Minv = np.linalg.solve(M_cat02, (M_cat02.T / D_RGB).T)
    else:
        # inverse as given in the standard
        M_cat02_inv = np.array(
            [
                [+1.096124, -0.278869, 0.182745],
                [+0.454369, +0.473533, 0.072098],
                [-0.009628, -0.005698, 1.015326],
            ]
        )
        M = M_cat02_inv @ (M_cat02.T * D_RGB).T
        Minv = M_cat02_inv @ (M_cat02.T / D_RGB).T

    return M, Minv
