import numpy as np

from ._cielab import CIELAB
from ._color_space import ColorSpace


class DIN99(ColorSpace):
    def __init__(self, k_E=1.0, k_CH=1.0):
        super().__init__("DIN99", ("L99", "a99", "b99"), 0)
        self.k_E = k_E
        self.k_CH = k_CH
        self.cielab = CIELAB()

    def from_xyz100(self, xyz):
        L, a, b = self.cielab.from_xyz100(xyz)
        L99 = 105.51 * np.log(1 + 0.0158 * L) / self.k_E

        sin16 = np.sin(np.radians(16))
        cos16 = np.cos(np.radians(16))

        e = a * cos16 + b * sin16
        f = 0.7 * (-a * sin16 + b * cos16)

        G = np.hypot(e, f)

        k = np.log(1 + 0.045 * G) / 0.045

        a99 = k * e / G
        b99 = k * f / G

        a99 = np.nan_to_num(a99, nan=0.0)
        b99 = np.nan_to_num(b99, nan=0.0)

        return np.array([L99, a99, b99])

    def to_xyz100(self, lab99):
        L99, a99, b99 = lab99
        C99 = np.hypot(a99, b99)
        G = (np.exp(0.045 * C99 * self.k_CH * self.k_E) - 1) / 0.045

        hp = np.hypot(a99, b99)
        cosarctan = a99 / hp
        sinarctan = b99 / hp
        cosarctan = np.nan_to_num(cosarctan, nan=0.0)
        sinarctan = np.nan_to_num(sinarctan, nan=0.0)

        e = G * cosarctan
        f = G * sinarctan

        sin16 = np.sin(np.radians(16))
        cos16 = np.cos(np.radians(16))
        a = e * cos16 - f / 0.7 * sin16
        b = e * sin16 + f / 0.7 * cos16

        L = (np.exp(L99 * self.k_E / 105.51) - 1) / 0.0158

        lab = np.array([L, a, b])
        return self.cielab.to_xyz100(lab)
