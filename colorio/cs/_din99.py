import numpy as np

from ._cielab import CIELAB
from ._color_space import ColorSpace


class DIN99(ColorSpace):
    def __init__(self, k_E=1.0, k_CH=1.0, variant=None):
        # variants from
        #
        # G. Cui, M.R. Luo, B. Rigg, G. Roesler, K. Witt,
        # Uniform colour spaces based on the DIN99 colour-difference formula
        # <https://doi.org/10.1002/col.10066>.
        name = "DIN99"
        if variant is not None:
            name += variant
        super().__init__(
            f"DIN99{variant}", (f"L99{variant}", f"a99{variant}", f"b99{variant}"), 0
        )
        self.k_E = k_E
        self.k_CH = k_CH
        self.cielab = CIELAB()

        if variant is None:
            self.p = [105.51, 0.0158, 16.0, 0.7, 200 / 9, 9 / 200, 0.0]
        elif variant == "b":
            self.p = [303.67, 0.0039, 26.0, 0.83, 23.0, 0.075, 26 / 180 * np.pi]
        elif variant == "c":
            self.p = [317.65, 0.0037, 0.0, 0.94, 23.0, 0.066, 0.0]
        else:
            assert variant == "d"
            self.p = [325.22, 0.0036, 50.0, 1.14, 22.5, 0.06, 50 / 180 * np.pi]

    def from_xyz100(self, xyz):
        L, a, b = self.cielab.from_xyz100(xyz)
        L99 = self.p[0] * np.log(1 + self.p[1] * L) / self.k_E

        sin = np.sin(np.radians(self.p[2]))
        cos = np.cos(np.radians(self.p[2]))

        e = a * cos + b * sin
        f = self.p[3] * (-a * sin + b * cos)

        G = np.hypot(e, f)

        C99 = self.p[4] * np.log(1 + self.p[5] * G)

        # h99 = np.arctan2(f, e) + self.p[6]
        # cos_h99 = np.cos(h99)
        # sin_h99 = np.sin(h99)

        cosarctan_f_e = e / G
        sinarctan_f_e = f / G
        sin_p6 = np.sin(self.p[6])
        cos_p6 = np.cos(self.p[6])
        cos_h99 = cosarctan_f_e * cos_p6 - sinarctan_f_e * sin_p6
        sin_h99 = sinarctan_f_e * cos_p6 + cosarctan_f_e * sin_p6

        a99 = C99 * cos_h99
        b99 = C99 * sin_h99

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
