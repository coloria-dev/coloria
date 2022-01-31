"""
https://en.wikipedia.org/wiki/Color_difference#CIE94
"""
import numpy as np
from numpy.typing import ArrayLike


# parameters:
# graphic arts:
#   k_L = 1
#   K_1 = 0.045
#   K_2 = 0.015
#
# textiles:
#   k_L = 2
#   K_1 = 0.048
#   K_2 = 0.014
def cie94(
    lab1: ArrayLike,
    lab2: ArrayLike,
    k_L: float = 1.0,
    K_1: float = 0.045,
    K_2: float = 0.015,
) -> np.ndarray:
    l1, a1, b1 = np.asarray(lab1)
    l2, a2, b2 = np.asarray(lab2)

    dL = l1 - l2
    C1 = np.sqrt(a1**2 + b1**2)
    C2 = np.sqrt(a2**2 + b2**2)
    dC = C1 - C2
    da = a1 - a2
    db = b1 - b2
    # dH2 is mathematically >=0, but round-off can lead to small negatives
    dH2 = da**2 + db**2 - dC**2

    S_L = 1.0
    S_C = 1 + K_1 * C1
    S_H = 1 + K_2 * C1

    k_C = 1.0
    k_H = 1.0
    dE = np.sqrt((dL / k_L / S_L) ** 2 + (dC / k_C / S_C) ** 2 + dH2 / (k_H * S_H) ** 2)
    return dE
