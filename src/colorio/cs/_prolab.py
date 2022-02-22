"""
Ivan A. Konovalenko, Anna A. Smagina, Dmitry P. Nikolaev, Petr P. Nikolaev,
ProLab: perceptually uniform projective colour coordinate system,
https://arxiv.org/abs/2012.07653
"""
import npx
import numpy as np
from numpy.typing import ArrayLike

from ..illuminants import whitepoints_cie1931
from ._color_space import ColorSpace
from ._helpers import register


class PROLAB(ColorSpace):
    name = "proLab"
    labels = ("L", "a", "b")
    k0 = 0

    def __init__(self, whitepoint=whitepoints_cie1931["D65"]):
        # The matrix Q in the article is
        #
        #  Q = ( Q 0 )
        #      ( q 1 )
        #
        self.Q = np.array(
            [
                [75.54, 486.66, 167.39],
                [617.72, -595.45, -22.27],
                [48.34, 194.94, -243.28],
            ]
        )
        self.q = np.array([0.7554, 3.8666, 1.6739])
        self.Qinv = np.linalg.inv(self.Q)
        self.whitepoint_xyz100 = np.asarray(whitepoint)
        self.whitepoint = np.array([100.0, 0.0, 0.0])

        # P is Q with whitepoint normalization
        # self.P = np.array([
        #     [79.4725, 486.6610, 153.7311],
        #     [649.9038, -595.4477, -20.4498],
        #     [50.8625, 194.9377, -223.4334],
        #     ])
        # self.p = np.array([0.7947, 3.8666, 1.5373])

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        xyz = np.asarray(xyz)
        xyz = (xyz.T / self.whitepoint_xyz100).T
        return npx.dot(self.Q, xyz) / (npx.dot(self.q, xyz) + 1)

    def to_xyz100(self, lab: ArrayLike) -> np.ndarray:
        y = npx.dot(self.Qinv, lab)
        xyz = y / (1 - npx.dot(self.q, y))
        xyz = (xyz.T * self.whitepoint_xyz100).T
        return xyz


register("prolab", PROLAB())
