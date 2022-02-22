"""
https://de.wikipedia.org/wiki/DIN99-Farbraum
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from ._cielab import CIELAB
from ._color_space import ColorSpace
from ._helpers import register


class DIN99(ColorSpace):
    name = "DIN99"
    labels = ("L99", "a99", "b99")
    k0 = 0

    def __init__(self, k_E: float = 1.0, k_CH: float = 1.0, variant: str | None = None):
        # variants from
        #
        # G. Cui, M.R. Luo, B. Rigg, G. Roesler, K. Witt,
        # Uniform colour spaces based on the DIN99 colour-difference formula
        # <https://doi.org/10.1002/col.10066>.
        self.k_E = k_E
        self.k_CH = k_CH
        self.cielab = CIELAB()

        if variant is None:
            self.p = [105.51, 0.0158, 16.0, 0.7, 200 / 9, 9 / 200, 0.0]
        elif variant == "b":
            self.p = [303.67, 0.0039, 26.0, 0.83, 23.0, 0.075, 26.0]
        elif variant == "c":
            self.p = [317.65, 0.0037, 0.0, 0.94, 23.0, 0.066, 0.0]
        else:
            assert variant == "d"
            self.p = [325.22, 0.0036, 50.0, 1.14, 22.5, 0.06, 50.0]

        self.sin_p2 = np.sin(np.radians(self.p[2]))
        self.cos_p2 = np.cos(np.radians(self.p[2]))
        self.sin_p6 = np.sin(np.radians(self.p[6]))
        self.cos_p6 = np.cos(np.radians(self.p[6]))

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        L, a, b = self.cielab.from_xyz100(xyz)
        L99 = self.p[0] * np.log(1 + self.p[1] * L) / self.k_E

        e = a * self.cos_p2 + b * self.sin_p2
        f = self.p[3] * (-a * self.sin_p2 + b * self.cos_p2)

        G = np.hypot(e, f)

        C99 = self.p[4] * np.log(1 + self.p[5] * G)

        # h99 = np.arctan2(f, e) + self.p[6]
        # cos_h99 = np.cos(h99)
        # sin_h99 = np.sin(h99)

        cosarctan_f_e = np.zeros_like(G)
        np.divide(e, G, out=cosarctan_f_e, where=G != 0.0)

        sinarctan_f_e = np.zeros_like(G)
        np.divide(f, G, out=sinarctan_f_e, where=G != 0.0)

        cos_h99 = self.cos_p6 * cosarctan_f_e - self.sin_p6 * sinarctan_f_e
        sin_h99 = self.sin_p6 * cosarctan_f_e + self.cos_p6 * sinarctan_f_e

        a99 = C99 * cos_h99
        b99 = C99 * sin_h99

        return np.array([L99, a99, b99])

    def to_xyz100(self, lab99: ArrayLike) -> np.ndarray:
        L99, a99, b99 = np.asarray(lab99)
        C99 = np.hypot(a99, b99)
        G = (np.exp(C99 / self.p[4] * self.k_CH * self.k_E) - 1) / self.p[5]

        cos_h99 = np.zeros_like(C99)
        np.divide(a99, C99, out=cos_h99, where=C99 != 0.0)
        sin_h99 = np.zeros_like(C99)
        np.divide(b99, C99, out=sin_h99, where=C99 != 0.0)

        cosarctan_f_e = self.cos_p6 * cos_h99 + self.sin_p6 * sin_h99
        sinarctan_f_e = -self.sin_p6 * cos_h99 + self.cos_p6 * sin_h99

        e = G * cosarctan_f_e
        f = G * sinarctan_f_e

        a = e * self.cos_p2 - f / self.p[3] * self.sin_p2
        b = e * self.sin_p2 + f / self.p[3] * self.cos_p2

        L = (np.exp(L99 * self.k_E / self.p[0]) - 1) / self.p[1]

        return self.cielab.to_xyz100([L, a, b])


register("din99", DIN99())
register("din99b", DIN99(variant="b"))
register("din99c", DIN99(variant="c"))
register("din99d", DIN99(variant="d"))
