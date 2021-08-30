"""
S. Bianco, R. Schettini,
Two New von Kries Based Chromatic Adaptation Transforms Found by Numerical Optimization,
2010,
<https://doi.org/10.1002/col.20573>
"""
import numpy as np
from numpy.typing import ArrayLike

from .von_kries import _von_kries_type


def bianco_schettini(whitepoint_source: ArrayLike, whitepoint_target: ArrayLike):
    M = np.array(
        [
            [0.8752, 0.2787, -0.1539],
            [-0.8904, 1.8709, 0.0195],
            [-0.0061, 0.0162, 0.9899],
        ]
    )
    return _von_kries_type(M, whitepoint_source, whitepoint_target)


def bianco_schettini_pos(whitepoint_source: ArrayLike, whitepoint_target: ArrayLike):
    M = np.array(
        [
            [0.6489, 0.3915, -0.0404],
            [-0.3775, 1.3055, 0.0720],
            [-0.0271, 0.0888, 0.9383],
        ]
    )
    return _von_kries_type(M, whitepoint_source, whitepoint_target)
