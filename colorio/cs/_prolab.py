# Ivan A. Konovalenko, Anna A. Smagina, Dmitry P. Nikolaev, Petr P. Nikolaev,
# ProLab: perceptually uniform projective colour coordinate system,
# https://arxiv.org/abs/2012.07653
import numpy as np

from .._linalg import dot
from ..illuminants import whitepoints_cie1931
from ._color_space import ColorSpace


class PROLAB(ColorSpace):
    def __init__(self, whitepoint=whitepoints_cie1931["D65"]):
        super().__init__("proLab", ("L", "a", "b"), 0)
        # The matrix Q in the article is
        #
        #  Q = ( Q 0 )
        #      ( q 1 )
        #

        # full-precision Q from
        # https://github.com/konovalenko-iitp/proLab/blob/main/matlab/proLab_param.m
        self.Q = 100 * np.array(
            [
                [0.755362206159585, 4.866609636500914, 1.673869831378788],
                [6.177140886995143, -5.954477384591136, -0.222663502404010],
                [0.483432883468812, 1.949376826435850, -2.432809709904663],
            ]
        )
        self.q = 100 * np.array(
            [0.007553621984060, 0.038666096468163, 0.016738698288170]
        )
        self.Qinv = np.linalg.inv(self.Q)
        self.wp = whitepoint

        # P is Q with whitepoint normalization
        # self.P = np.array([
        #     [79.4725, 486.6610, 153.7311],
        #     [649.9038, -595.4477, -20.4498],
        #     [50.8625, 194.9377, -223.4334],
        #     ])
        # self.p = np.array([0.7947, 3.8666, 1.5373])

    def from_xyz100(self, xyz):
        xyz = np.asarray(xyz)
        xyz = (xyz.T / self.wp).T
        return dot(self.Q, xyz) / (dot(self.q, xyz) + 1)

    def to_xyz100(self, lab):
        y = dot(self.Qinv, lab)
        xyz = y / (1 - dot(self.q, y))
        xyz = (xyz.T * self.wp).T
        return xyz
