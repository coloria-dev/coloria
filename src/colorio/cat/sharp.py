"""
See, e.g., <http://dx.doi.org/10.1002/col.20573>.
"""
import numpy as np
from numpy.typing import ArrayLike


def sharp(whitepoint_source: ArrayLike, whitepoint_target: ArrayLike):
    M = np.array(
        [
            [1.2694, 0.0988, -0.1706],
            [-0.8364, 1.8006, 0.0357],
            [0.0297, -0.0315, 1.0018],
        ]
    )

    whitepoint_source = np.asarray(whitepoint_source)
    whitepoint_target = np.asarray(whitepoint_target)

    rgb_rw = M @ whitepoint_target
    rgb_w = M @ whitepoint_source
    abc = rgb_rw / rgb_w

    A = np.linalg.solve(M, (M.T * abc).T)
    Ainv = np.linalg.solve(M, (M.T / abc).T)
    return A, Ainv
