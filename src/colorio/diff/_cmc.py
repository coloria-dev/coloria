"""
https://en.wikipedia.org/wiki/Color_difference#CMC_l:c_(1984)
"""
import numpy as np
from numpy.typing import ArrayLike

from ..cs import ColorCoordinates, convert


# Acceptability: l=2, c=1 Perceptability: l=1, c=1
def cmc(lab1: ArrayLike, lab2: ArrayLike, l: float = 2.0, c: float = 1.0) -> np.ndarray:
    lab1 = ColorCoordinates(lab1, "cielab")
    lab2 = ColorCoordinates(lab2, "cielab")

    lch1 = convert(lab1, "cielch")
    lch2 = convert(lab2, "cielch")

    L1, C1, h1 = lch1.data
    L2, C2, _ = lch2.data

    F = np.sqrt(C1**4 / (C1**4 + 1900))

    idx = (164 <= h1) & (h1 <= 345)
    #
    p1 = np.empty_like(h1)
    p1[idx] = 0.56
    p1[~idx] = 0.36
    #
    p2 = np.empty_like(h1)
    p2[idx] = 0.2
    p2[~idx] = 0.4
    #
    offset = np.empty_like(h1)
    offset[idx] = 168
    offset[~idx] = 35

    T = p1 + np.abs(p2 * np.cos(np.radians(h1 + offset)))

    S_L = np.empty_like(L1)
    idx = L1 < 16
    S_L[idx] = 0.511
    S_L[~idx] = (0.040975 * L1[~idx]) / (1 + 0.01765 * L1[~idx])

    S_C = 0.0638 * C1 / (1 + 0.0131 * C1) + 0.638
    S_H = S_C * (F * T + 1 - F)

    _, a1, b1 = lab1.data
    _, a2, b2 = lab2.data
    dC = C1 - C2
    da = a1 - a2
    db = b1 - b2
    dHab2 = da**2 + db**2 - dC**2

    dE = np.sqrt(((L2 - L1) / l / S_L) ** 2 + (dC / c / S_C) ** 2 + dHab2 / S_H**2)
    return dE
