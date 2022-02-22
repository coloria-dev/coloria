import npx
import numpy as np
from numpy.typing import ArrayLike

from .._exceptions import ColorioError
from ._color_space import ColorSpace
from ._helpers import register


class JzAzBz(ColorSpace):
    """
    Muhammad Safdar, Guihua Cui, Youn Jin Kim, and Ming Ronnier Luo,
    Perceptually uniform color space for image signals including high dynamic range and
    wide gamut,
    Optics Express Vol. 25, Issue 13, pp. 15131-15151 (2017),
    <https://doi.org/10.1364/OE.25.015131>.
    """

    name = "J_z a_z b_z"
    labels = ("J_z", "a_z", "b_z")
    k0 = 0

    def __init__(self):
        self.b = 1.15
        self.g = 0.66
        self.c1 = 3424 / 2**12
        self.c2 = 2413 / 2**7
        self.c3 = 2392 / 2**7
        self.n = 2610 / 2**14
        self.p = 1.7 * 2523 / 2**5
        self.d = -0.56
        self.d0 = 1.6295499532821566e-11

        self.M1 = np.array(
            [
                [0.41478972, 0.579999, 0.0146480],
                [-0.2015100, 1.120649, 0.0531008],
                [-0.0166008, 0.264800, 0.6684799],
            ]
        )
        self.M2 = np.array(
            [
                [0.5, 0.5, 0],
                [3.524000, -4.066708, +0.542708],
                [0.199076, +1.096799, -1.295875],
            ]
        )

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        # x, y, z = (xyz.T / self.whitepoint).T
        # In nit units, ranging from 0 to 10000?
        x, y, z = np.asarray(xyz)
        x_ = self.b * x - (self.b - 1) * z
        y_ = self.g * y - (self.g - 1) * x
        lms = npx.dot(self.M1, [x_, y_, z])
        lms_ = (
            (self.c1 + self.c2 * (lms / 10000) ** self.n)
            / (1 + self.c3 * (lms / 10000) ** self.n)
        ) ** self.p
        iz, az, bz = npx.dot(self.M2, lms_)
        jz = (1 + self.d) * iz / (1 + self.d * iz) - self.d0
        return np.array([jz, az, bz])

    def to_xyz100(self, jzazbz: ArrayLike) -> np.ndarray:
        jz, az, bz = np.asarray(jzazbz)
        iz = (jz + self.d0) / (1 + self.d - self.d * (jz + self.d0))
        lms_ = npx.solve(self.M2, np.array([iz, az, bz]))
        if np.any(lms_ < 0.0):
            raise ColorioError("Illegal LMS value.")

        lms = 10000 * (
            (self.c1 - lms_ ** (1 / self.p))
            / (self.c3 * lms_ ** (1 / self.p) - self.c2)
        ) ** (1 / self.n)
        x_, y_, z_ = npx.solve(self.M1, lms)
        x = (x_ + (self.b - 1) * z_) / self.b
        y = (y_ + (self.g - 1) * x) / self.g
        # return (np.array([x, y, z_]).T * self.whitepoint).T
        return np.array([x, y, z_])


register("jzazbz", JzAzBz())
