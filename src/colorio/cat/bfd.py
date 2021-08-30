"""
Bradford Chromatic Adaptation Transform (BFD CAT)
"""
import numpy as np
from numpy.typing import ArrayLike


def bfd(whitepoint_source: ArrayLike, whitepoint_target: ArrayLike):
    M = np.array(
        [
            [+0.8951, +0.2664, -0.1614],
            [-0.7502, +1.7135, +0.0367],
            [+0.0389, -0.0685, +1.0296],
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
