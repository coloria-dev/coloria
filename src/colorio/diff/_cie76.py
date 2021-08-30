"""
https://en.wikipedia.org/wiki/Color_difference#CIE76
"""
import numpy as np
from numpy.typing import ArrayLike


def cie76(lab1: ArrayLike, lab2: ArrayLike) -> np.ndarray:
    diff = np.asarray(lab1) - np.asarray(lab2)
    return np.sqrt(np.einsum("i...,i...->...", diff, diff))
