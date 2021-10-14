"""
https://en.wikipedia.org/wiki/Color_difference#CIEDE2000
"""
import numpy as np
from numpy.typing import ArrayLike


def ciede2000(
    lab1: ArrayLike,
    lab2: ArrayLike,
    k_L: float = 1.0,
    k_C: float = 1.0,
    k_H: float = 1.0,
) -> np.ndarray:
    lab1 = np.asarray(lab1)
    lab2 = np.asarray(lab2)

    L1, a1, b1 = lab1
    L2, a2, b2 = lab2

    C1 = np.sqrt(a1 ** 2 + b1 ** 2)
    C2 = np.sqrt(a2 ** 2 + b2 ** 2)
    C_mean = (C1 + C2) / 2

    G = 0.5 * (1 - np.sqrt(C_mean ** 7 / (C_mean ** 7 + 25 ** 7)))
    a1p = (1 + G) * a1
    a2p = (1 + G) * a2

    C1p = np.sqrt(a1p ** 2 + b1 ** 2)
    C2p = np.sqrt(a2p ** 2 + b2 ** 2)

    h1p = np.degrees(np.arctan2(b1, a1p)) % 360
    h2p = np.degrees(np.arctan2(b2, a2p)) % 360
    # make sure dhp is in the (-180, 180) range
    hp_diff = h2p - h1p
    dhp = ((hp_diff + 180) % 360) - 180

    dLp = L2 - L1
    dCp = C2p - C1p
    dHp = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dhp / 2))

    Lp_mean = (L1 + L2) / 2
    Cp_mean = (C1p + C2p) / 2

    hp_mean = np.asarray((h1p + h2p) / 2)

    idx = np.abs(hp_diff) > 180
    hp_mean[idx] = (hp_mean[idx] - 180) % 360

    T = (
        1.0
        - 0.17 * np.cos(np.radians(hp_mean - 30))
        + 0.24 * np.cos(np.radians(2 * hp_mean))
        + 0.32 * np.cos(np.radians(3 * hp_mean + 6))
        - 0.20 * np.cos(np.radians(4 * hp_mean - 63))
    )
    dtheta = 30 * np.exp(-(((hp_mean - 275) / 25) ** 2))
    R_C = 2 * np.sqrt(Cp_mean ** 7 / (Cp_mean ** 7 + 25 ** 7))
    S_L = 1 + 0.015 * (Lp_mean - 50) ** 2 / np.sqrt(20 + (Lp_mean - 50) ** 2)
    S_C = 1 + 0.045 * Cp_mean
    S_H = 1 + 0.015 * Cp_mean * T
    R_T = -np.sin(np.radians(2 * dtheta)) * R_C

    dE00 = np.sqrt(
        (dLp / k_L / S_L) ** 2
        + (dCp / k_C / S_C) ** 2
        + (dHp / k_H / S_H) ** 2
        + R_T * (dCp / k_C / S_C) * (dHp / k_H / S_H)
    )
    return dE00
