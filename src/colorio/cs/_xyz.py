import numpy as np
from numpy.typing import ArrayLike

from .._exceptions import ColorioError
from ._color_space import ColorSpace


class XYZ(ColorSpace):
    """XYZ colorspace."""

    def __init__(self, scaling: float):
        if scaling not in [1, 100]:
            raise ColorioError("scaling needs to be 1 or 100.")

        self.scaling = scaling
        super().__init__(f"XYZ_{self.scaling}", ("X", "Y", "Z"), None)

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        xyz = np.asarray(xyz)
        if self.scaling == 100:
            return xyz

        return xyz / 100

    def to_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        xyz = np.asarray(xyz)
        if self.scaling == 100:
            return xyz
        return xyz * 100


class XYZ1(XYZ):
    def __init__(self):
        super().__init__(1)


class XYZ100(XYZ):
    def __init__(self):
        super().__init__(100)
