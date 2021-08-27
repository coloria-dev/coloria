import numpy as np
from numpy.typing import ArrayLike

from ..illuminants import whitepoints_cie1931
from ._cielab import f, finv
from ._color_space import ColorSpace


class CIELUV(ColorSpace):
    name = "CIELUV"
    labels = ("L*", "u*", "v*")
    k0 = 0
    is_origin_well_defined = False

    def __init__(self, whitepoint: ArrayLike = whitepoints_cie1931["D65"]):
        self.whitepoint_xyz100 = np.asarray(whitepoint)
        self.whitepoint = np.array([100.0, 0.0, 0.0])

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        xyz = np.asarray(xyz)
        L = f(xyz[1] / self.whitepoint_xyz100[1])

        x, y, z = xyz
        p = x + 15 * y + 3 * z
        u = 4 * x / p
        v = 9 * y / p

        wx, wy, wz = self.whitepoint_xyz100
        q = wx + 15 * wy + 3 * wz
        un = 4 * wx / q
        vn = 9 * wy / q
        return np.array([L, 13 * L * (u - un), 13 * L * (v - vn)])

    def to_xyz100(self, luv: ArrayLike) -> np.ndarray:
        L, u, v = np.asarray(luv)

        wx, wy, wz = self.whitepoint_xyz100
        q = wx + 15 * wy + 3 * wz
        un = 4 * wx / q
        vn = 9 * wy / q

        uu = u / (13 * L) + un
        vv = v / (13 * L) + vn

        Y = wy * finv(L)
        X = Y * 9 * uu / (4 * vv)
        Z = Y * (12 - 3 * uu - 20 * vv) / (4 * vv)
        return np.array([X, Y, Z])
