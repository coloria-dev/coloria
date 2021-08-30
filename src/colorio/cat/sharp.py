"""
See, e.g., <http://dx.doi.org/10.1002/col.20573>.
"""
import numpy as np
from numpy.typing import ArrayLike

from .von_kries import _von_kries_type


def sharp(whitepoint_source: ArrayLike, whitepoint_target: ArrayLike):
    M = np.array(
        [
            [1.2694, 0.0988, -0.1706],
            [-0.8364, 1.8006, 0.0357],
            [0.0297, -0.0315, 1.0018],
        ]
    )
    return _von_kries_type(M, whitepoint_source, whitepoint_target)
