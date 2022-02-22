import npx
import numpy as np
from numpy.typing import ArrayLike

from ._color_space import ColorSpace
from ._hdr import HdrLinear
from ._helpers import register


class ICtCp(ColorSpace):
    """
    ICtCp color model.
    <https://en.wikipedia.org/wiki/ICtCp>
    """

    name = "IC_TC_P"
    labels = ("I", "C_T", "C_P")
    k0 = 0

    def __init__(self):
        self.M1 = (
            np.array([[1688, 2146, 262], [683, 2951, 462], [99, 309, 3688]]) / 4096
        )

        # From <https://doi.org/10.5594/SMPTE.ST2084.2014>
        self.m1 = 2610 / 4096 / 4
        self.m2 = 2523 / 4096 * 128
        self.c1 = 3424 / 4096
        self.c2 = 2413 / 4096 * 32
        self.c3 = 2392 / 4096 * 32

        self.M2 = (
            np.array([[2048, 2048, 0], [6610, -13613, 7003], [17933, -17390, -543]])
            / 4096
        )

        self._hdr = HdrLinear()

    def from_rec2100(self, rgb: ArrayLike) -> ArrayLike:
        lms = npx.dot(self.M1, rgb)

        lms_ = (
            (self.c1 + self.c2 * lms**self.m1) / (1 + self.c3 * lms**self.m1)
        ) ** self.m2

        ictcp = npx.dot(self.M2, lms_)
        return ictcp

    def to_rec2100(self, ictcp: ArrayLike) -> ArrayLike:
        lms_ = npx.solve(self.M2, ictcp)

        t = lms_ ** (1 / self.m2) - self.c1
        # This next line is part of the model, but really it shouldn't occur for sane
        # input data.
        # t[t < 0] = 0.0
        lms = (t / (self.c2 - self.c3 * lms_ ** (1 / self.m2))) ** (1 / self.m1)

        rgb = npx.solve(self.M1, lms)
        return rgb

    def from_xyz100(self, xyz100: ArrayLike) -> ArrayLike:
        return self.from_rec2100(self._hdr.from_xyz100(xyz100))

    def to_xyz100(self, ictcp: ArrayLike) -> ArrayLike:
        return self._hdr.to_xyz100(self.to_rec2100(ictcp))


register("ictcp", ICtCp())
