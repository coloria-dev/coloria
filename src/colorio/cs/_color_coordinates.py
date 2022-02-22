from __future__ import annotations

from copy import deepcopy

import numpy as np
from numpy.typing import ArrayLike

from ._color_space import ColorSpace
from ._helpers import string_to_cs


class ColorCoordinates:
    def __init__(self, data: ArrayLike, color_space: ColorSpace | str):
        self.data = np.asarray(data)
        if self.data.shape[0] != 3:
            raise ValueError(f"Input data needs shape [3,...], got {self.data.shape}.")

        if isinstance(color_space, str):
            color_space = string_to_cs(color_space)

        self.color_space = color_space

    def __repr__(self) -> str:
        return (
            "<ColorCoordinates, "
            + f"{self.color_space.name}, "
            + f"data.shape={self.data.shape}>"
        )

    def __mul__(self, alpha):
        return ColorCoordinates(alpha * self.data, self.color_space)

    __rmul__ = __mul__

    def __add__(self, other):
        if isinstance(other, ColorCoordinates):
            if self.color_space.name != other.color_space.name:
                raise ValueError(
                    "Color spaces not equal "
                    f"({self.color_space.name} != {other.color_space.name})"
                )
            return ColorCoordinates(self.data + other.data, self.color_space)

        # fallback for int/float etc
        return ColorCoordinates(self.data + other, self.color_space)

    __radd__ = __add__

    def __eq__(self, other):
        if isinstance(other, ColorCoordinates):
            if self.color_space.name != other.color_space.name:
                raise ValueError(
                    "Color spaces not equal "
                    f"({self.color_space.name} != {other.color_space.name})"
                )
            return self.data == other.data

        # fallback for int/float etc
        return self.data == other

    def __lt__(self, other):
        return self.data < other

    def __le__(self, other):
        return self.data <= other

    def __gt__(self, other):
        return self.data > other

    def __ge__(self, other):
        return self.data >= other

    def copy(self):
        return deepcopy(self)

    def convert(self, cs: ColorSpace | str, **kwargs) -> None:
        if isinstance(cs, str):
            cs = string_to_cs(cs)

        if cs == self.color_space:
            return

        self.data = cs.from_xyz100(self.color_space.to_xyz100(self.data), **kwargs)
        self.color_space = cs

    @property
    def lightness(self) -> np.ndarray:
        assert self.color_space.k0 is not None
        return self.data[self.color_space.k0]

    @property
    def hue(self) -> np.ndarray:
        assert self.color_space.k0 is not None
        hue_idx = np.array([True, True, True])
        hue_idx[self.color_space.k0] = False
        return self.data[hue_idx]


def convert(
    coords: ColorCoordinates, cs: ColorSpace | str, **kwargs
) -> ColorCoordinates:
    out = coords.copy()
    out.convert(cs, **kwargs)
    return out
