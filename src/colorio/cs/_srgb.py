"""
https://en.wikipedia.org/wiki/SRGB
"""
from __future__ import annotations

import npx
import numpy as np
from numpy.typing import ArrayLike

from ..illuminants import whitepoints_cie1931
from ._color_space import ColorSpace
from ._helpers import register


def _xyy_to_xyz100(xyy):
    x, y, Y = xyy
    return np.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100


class SRGBlinear(ColorSpace):
    """Rec. 709 SRGB."""

    def __init__(self, default_mode: str = "error", whitepoint_correction: bool = True):
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

        self.default_mode = default_mode
        self.name = "sRGB (linear)"

        # np.linalg.inv(self.invM) is the matrix in the spec:
        # M = np.array([
        #     [+3.2406255, -1.537208, -0.4986286],
        #     [-0.9689307, +1.8757561, +0.0415175],
        #     [+0.0557101, -0.2040211, +1.0569959],
        # ])
        # self.invM = np.linalg.inv(M)
        self.labels = ["R", "G", "B"]

    def from_xyz100(self, xyz: ArrayLike, mode: str | None = None) -> np.ndarray:
        if mode is None:
            mode = self.default_mode
        # https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
        # https://www.color.org/srgb.pdf
        out = npx.solve(self.invM, xyz) / 100

        if mode == "error":
            if np.any(out < 0) or np.any(out > 1):
                raise ValueError(
                    "Not all XYZ values could be converted to legal sRGB. "
                    + 'Try with `mode="clip"`.'
                )
        elif mode == "ignore":
            pass
        elif mode == "nan":
            out[out < 0] = np.nan
            out[out > 1] = np.nan
        else:
            assert mode == "clip"
            out = out.clip(0.0, 1.0)

        return out

    def to_xyz100(self, srgb1_linear: ArrayLike) -> np.ndarray:
        # Note: The Y value is often used for grayscale conversion.
        # 0.2126 * R_linear + 0.7152 * G_linear + 0.0722 * B_linear
        return 100 * npx.dot(self.invM, srgb1_linear)


class SRGB1(ColorSpace):
    def __init__(self, default_mode: str = "error"):
        self._srgb_linear = SRGBlinear(default_mode=default_mode)
        self.name = "sRGB-1"

    def from_xyz100(self, xyz: ArrayLike, mode: str | None = None) -> np.ndarray:
        srgb = self._srgb_linear.from_xyz100(xyz, mode=mode)

        a = 0.055
        is_smaller = srgb <= 0.0031308
        srgb[is_smaller] *= 12.92
        srgb[~is_smaller] = (1 + a) * srgb[~is_smaller] ** (1 / 2.4) - a
        return srgb

    def to_xyz100(self, coords: ArrayLike) -> np.ndarray:
        coords = np.asarray(coords)

        a = 0.055
        # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
        is_smaller = coords <= 0.040449936  # 12.92 * 0.0031308
        coords[is_smaller] /= 12.92
        coords[~is_smaller] = ((coords[~is_smaller] + a) / (1 + a)) ** 2.4

        return self._srgb_linear.to_xyz100(coords)


class SRGB255(ColorSpace):
    def __init__(self, default_mode: str = "error"):
        self._srgb1 = SRGB1(default_mode=default_mode)
        self.name = "sRGB-255"

    def from_xyz100(self, xyz: ArrayLike, mode: str | None = None) -> np.ndarray:
        return 255 * self._srgb1.from_xyz100(xyz, mode=mode)

    def to_xyz100(self, coords: ArrayLike) -> np.ndarray:
        return self._srgb1.to_xyz100(np.asarray(coords) / 255)


class SRGBhex(ColorSpace):
    def __init__(self, default_mode: str = "error", prepend: str = "#"):
        self._srgb255 = SRGB255(default_mode=default_mode)
        self.name = "sRGB-hex"
        self.prepend = prepend

    def from_xyz100(self, xyz: ArrayLike, mode: str | None = None) -> np.ndarray:
        rgb255 = self._srgb255.from_xyz100(xyz, mode=mode)

        # round to closest int
        rgb255_rounded = np.around(rgb255).astype(int)

        if mode == "error":
            raise ValueError(
                "Rounding in sRGB-hex conversion "
                + f"from\n\n{rgb255.tolist()}\nto\n{rgb255_rounded.tolist()}\n"
            )

        # convert to hex, preserve shape
        shape = rgb255_rounded.shape
        assert shape[0] == 3
        rgb255_rounded = rgb255_rounded.reshape(3, -1)
        hex_vals = np.array(
            [self.prepend + f"{r:02x}{g:02x}{b:02x}" for r, g, b in rgb255_rounded.T]
        ).reshape(shape[1:])

        return hex_vals

    def to_xyz100(self, coords: ArrayLike) -> np.ndarray:
        def _string_to_rgb255(string):
            return [
                int(string[0:2], 16),
                int(string[2:4], 16),
                int(string[4:6], 16),
            ]

        coords = np.asarray(coords)

        shape = coords.shape

        srgb255 = [
            _string_to_rgb255(coord.item()[len(self.prepend) :])
            for coord in coords.reshape(-1)
        ]
        return self._srgb255.to_xyz100(np.asarray(srgb255).T).reshape(3, *shape)


register("srgblinear", SRGBlinear())
register("srgb1", SRGB1())
register("srgb255", SRGB255())
register("srgbhex", SRGBhex())
