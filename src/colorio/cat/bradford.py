"""
Bradford Chromatic Adaptation Transform (BFD CAT)
"""
import numpy as np
from numpy.typing import ArrayLike

from .von_kries import _von_kries_type


def bradford(whitepoint_source: ArrayLike, whitepoint_target: ArrayLike):
    M = np.array(
        [
            [+0.8951, +0.2664, -0.1614],
            [-0.7502, +1.7135, +0.0367],
            [+0.0389, -0.0685, +1.0296],
        ]
    )
    return _von_kries_type(M, whitepoint_source, whitepoint_target)
