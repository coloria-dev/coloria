import npx
import numpy as np
from numpy.typing import ArrayLike

from ..illuminants import whitepoints_cie1931


def _xyy_to_xyz100(xyy):
    x, y, Y = xyy
    return np.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100


class SrgbLinear:
    """Rec. 709 SRGB."""

    def __init__(self, whitepoint_correction: bool = True):
        # The standard actually gives the values in terms of M, but really inv(M) is a
        # direct derivative of the primary specification at
        # <https://en.wikipedia.org/wiki/SRGB>.
        primaries_xyy = np.array(
            [[0.64, 0.33, 0.2126], [0.30, 0.60, 0.7152], [0.15, 0.06, 0.0722]]
        )
        self.invM = _xyy_to_xyz100(primaries_xyy.T)

        if whitepoint_correction:
            # The above values are given only approximately, resulting in the fact that
            # SRGB(1.0, 1.0, 1.0) is only approximately mapped into the reference
            # whitepoint D65. Add a correction here.
            correction = whitepoints_cie1931["D65"] / np.sum(self.invM, axis=1)
            self.invM = (self.invM.T * correction).T

        self.invM /= 100

        # np.linalg.inv(self.invM) is the matrix in the spec:
        # M = np.array([
        #     [+3.2406255, -1.537208, -0.4986286],
        #     [-0.9689307, +1.8757561, +0.0415175],
        #     [+0.0557101, -0.2040211, +1.0569959],
        # ])
        # self.invM = np.linalg.inv(M)
        self.labels = ["R", "G", "B"]

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        # https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
        # https://www.color.org/srgb.pdf
        # TODO NaN the values smaller than 0 and larger than 1
        return npx.solve(self.invM, xyz) / 100

    def to_xyz100(self, srgb1_linear: ArrayLike) -> np.ndarray:
        # Note: The Y value is often used for grayscale conversion.
        # 0.2126 * R_linear + 0.7152 * G_linear + 0.0722 * B_linear
        return 100 * npx.dot(self.invM, srgb1_linear)

    def from_rgb1(self, srgb1: ArrayLike) -> np.ndarray:
        srgb_linear = np.array(srgb1, dtype=float)

        a = 0.055
        # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
        is_smaller = srgb_linear <= 0.040449936  # 12.92 * 0.0031308

        srgb_linear[is_smaller] /= 12.92
        srgb_linear[~is_smaller] = ((srgb_linear[~is_smaller] + a) / (1 + a)) ** 2.4
        return srgb_linear

    def to_rgb1(self, srgb_linear: ArrayLike) -> np.ndarray:
        a = 0.055
        srgb = np.copy(srgb_linear)
        is_smaller = srgb <= 0.0031308
        srgb[is_smaller] *= 12.92
        srgb[~is_smaller] = (1 + a) * srgb[~is_smaller] ** (1 / 2.4) - a
        return srgb

    def from_rgb255(self, srgb255: ArrayLike) -> np.ndarray:
        return self.from_rgb1(np.asarray(srgb255) / 255)

    def to_rgb255(self, srgb_linear: ArrayLike) -> np.ndarray:
        return 255 * self.to_rgb1(srgb_linear)
