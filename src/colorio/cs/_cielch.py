import numpy as np
from numpy.typing import ArrayLike

from ..illuminants import whitepoints_cie1931
from ._cielab import CIELAB
from ._color_space import ColorSpace


class CIELCH(ColorSpace):
    def __init__(self, whitepoint: ArrayLike = whitepoints_cie1931["D65"]):
        super().__init__("CIELCH", ("L", "C", "h"), 0)
        self.cielab = CIELAB(whitepoint=whitepoint)

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        L, a, b = self.cielab.from_xyz100(xyz)
        C = np.hypot(a, b)
        h = np.degrees(np.arctan2(b, a)) % 360
        return np.array([L, C, h])

    def to_xyz100(self, lch: ArrayLike) -> np.ndarray:
        L, C, h = np.asarray(lch)
        h_ = np.radians(h)
        lab = np.array([L, C * np.cos(h_), C * np.sin(h_)])
        return self.cielab.to_xyz100(lab)
