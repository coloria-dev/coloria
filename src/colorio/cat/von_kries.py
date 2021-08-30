"""
See

Ming Ronnier Luo,
CIE Chromatic Adaptation; Comparison of von Kries, CIELAB, CMCCAT97 and CAT02,
<https://doi.org/10.1007/978-3-642-27851-8_321-1>.
"""
import numpy as np
from numpy.typing import ArrayLike


def von_kries(
    whitepoint_source: ArrayLike,
    whitepoint_target: ArrayLike,
    exact_inversion: bool = True,
):
    judd = np.array(
        [
            [0.0, 1.0, 0.0],
            [-0.46, 1.36, 0.1],
            [0.0, 0.0, 1.0],
        ]
    )

    whitepoint_source = np.asarray(whitepoint_source)
    whitepoint_target = np.asarray(whitepoint_target)

    rgb_rw = judd @ whitepoint_target
    rgb_w = judd @ whitepoint_source
    abc = rgb_rw / rgb_w

    if exact_inversion:
        M = np.linalg.solve(judd, (judd.T * abc).T)
        Minv = np.linalg.solve(judd, (judd.T / abc).T)
    else:
        judd_inv_approx = np.array(
            [
                [2.954, -2.174, 0.22],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        M = judd_inv_approx @ (judd.T * abc).T
        Minv = judd_inv_approx @ (judd.T / abc).T

    return M, Minv
