from typing import Optional, Tuple

import numpy as np

from ._srgb import SrgbLinear


class ColorSpace:
    def __init__(
        self,
        name: str,
        labels: Tuple[str, str, str],
        k0: Optional[int],
        is_origin_well_defined: bool = True,
    ):
        self.name = name
        self.labels = np.asarray(labels)
        self.k0 = k0  # the index that corresponds to luminosity
        self.is_origin_well_defined = is_origin_well_defined
        self.srgb_linear = SrgbLinear()

    def __repr__(self):
        return f"<colorio color space {self.name}>"

    def to_xyz100(self, cs_coords):
        raise NotImplementedError("ColorSpace needs to implement to_xyz100()")

    def from_xyz100(self, cs_coords):
        raise NotImplementedError("ColorSpace needs to implement from_xyz100()")

    def to_rgb_linear(self, cs_coords):
        return self.srgb_linear.from_xyz100(self.to_xyz100(cs_coords))

    def to_rgb1(self, cs_coords):
        return self.srgb_linear.to_rgb1(self.to_rgb_linear(cs_coords))

    def to_rgb255(self, cs_coords):
        return self.srgb_linear.to_rgb255(self.to_rgb_linear(cs_coords))

    def from_rgb_linear(self, rgb_lin):
        return self.from_xyz100(self.srgb_linear.to_xyz100(rgb_lin))

    def from_rgb1(self, rgb1):
        return self.from_rgb_linear(self.srgb_linear.from_rgb1(rgb1))

    def from_rgb255(self, rgb255):
        return self.from_rgb_linear(self.srgb_linear.from_rgb255(rgb255))
