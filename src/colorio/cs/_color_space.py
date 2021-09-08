import numpy as np
from numpy.typing import ArrayLike

from ._srgb import SrgbLinear


class ColorSpace:
    is_origin_well_defined = True
    k0 = None
    srgb_linear = SrgbLinear()

    def __repr__(self):
        return f"<colorio color space {self.name}>"

    def to_xyz100(self, _):
        raise NotImplementedError("ColorSpace needs to implement to_xyz100()")

    def from_xyz100(self, _):
        raise NotImplementedError("ColorSpace needs to implement from_xyz100()")

    def to_rgb_hex(
        self, cs_coords: ArrayLike, mode: str = "error", prepend: str = "#"
    ) -> np.ndarray:
        xyz100 = self.to_xyz100(cs_coords)
        rgb_linear = self.srgb_linear.from_xyz100(xyz100, mode)
        return self.srgb_linear.to_rgb_hex(rgb_linear, prepend=prepend)
