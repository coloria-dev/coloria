import numpy as np
from numpy.typing import ArrayLike

from ._color_space import ColorSpace
from ._srgb import SrgbLinear


class ColorCoordinates:
    def __init__(self, data: ArrayLike, color_space: ColorSpace):
        self.data = np.asarray(data)
        if self.data.shape[0] != 3:
            raise ValueError(f"Input data needs shape [3,...], got {self.data.shape}.")
        self.color_space = color_space

    def __repr__(self):
        return (
            "<ColorCoordinates, "
            + f"{self.color_space.name}, "
            + f"data.shape={self.data.shape}>"
        )

    def convert(self, cs: ColorSpace):
        return ColorCoordinates(
            cs.from_xyz100(self.color_space.to_xyz100(self.data)), cs
        )

    @property
    def data_lightness(self):
        assert self.color_space.k0 is not None
        return self.data[self.color_space.k0]

    @property
    def data_hue(self):
        assert self.color_space.k0 is not None
        hue_idx = np.array([True, True, True])
        hue_idx[self.color_space.k0] = False
        return self.data[hue_idx]

    def get_rgb1(self, mode: str = "error"):
        s = SrgbLinear()
        return s.to_rgb1(
            s.from_xyz100(self.color_space.to_xyz100(self.data), mode=mode)
        )

    def get_rgb_hex(self, mode: str = "error", prepend: str = "#") -> np.ndarray:
        s = SrgbLinear()
        xyz100 = self.color_space.to_xyz100(self.data)
        rgb_linear = s.from_xyz100(xyz100, mode)
        return s.to_rgb_hex(rgb_linear, prepend=prepend)


def convert(coords: ColorCoordinates, cs: ColorSpace):
    return ColorCoordinates(
        cs.from_xyz100(coords.color_space.to_xyz100(coords.data)), cs
    )
