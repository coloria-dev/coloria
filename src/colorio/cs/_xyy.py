"""
xyY colorspace.

https://en.wikipedia.org/wiki/CIE_1931_color_space#CIE_xy_chromaticity_diagram_and_the_CIE_xyY_color_space
"""
import numpy as np
from numpy.typing import ArrayLike

from .._exceptions import ColorioError
from ._color_space import ColorSpace
from ._helpers import register


class XYY(ColorSpace):
    name = "xyY"
    labels = ("x", "y", "Y")
    k0 = 2
    is_origin_well_defined = False

    def __init__(self, Y_scaling: int):
        if Y_scaling not in [1, 100]:
            raise ColorioError("Y_scaling needs to be 1 or 100.")

        self.Y_scaling = Y_scaling

    def from_xyz100(self, xyz100: ArrayLike) -> np.ndarray:
        xyz100 = np.asarray(xyz100)
        if np.any(xyz100 < 0):
            raise ColorioError("Negative XYZ100 value.")

        xyz = xyz100
        if self.Y_scaling == 1:
            xyz = xyz100 / 100

        sum_xyz = np.sum(xyz, axis=0)
        x = xyz[0]
        y = xyz[1]
        return np.array([x / sum_xyz, y / sum_xyz, y])

    def to_xyz100(self, xyy: ArrayLike) -> np.ndarray:
        xyy = np.asarray(xyy)
        if np.any(xyy < 0):
            raise ColorioError("Negative xyY value.")
        x, y, Y = xyy
        out = np.array([Y / y * x, Y, Y / y * (1 - (x + y))])
        if self.Y_scaling == 1:
            out *= 100
        return out


class XYY1(XYY):
    name = "xyY1"
    labels = ("x", "y", "Y1")

    def __init__(self):
        super().__init__(1)


class XYY100(XYY):
    name = "xyY100"
    labels = ("x", "y", "Y100")

    def __init__(self):
        super().__init__(100)


register("xyy1", XYY1())
register("xyy100", XYY100())
