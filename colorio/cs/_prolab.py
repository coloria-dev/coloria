# Ivan A. Konovalenko, Anna A. Smagina, Dmitry P. Nikolaev, Petr P. Nikolaev,
# ProLab: perceptually uniform projective colour coordinate system,
# https://arxiv.org/abs/2012.07653
import numpy as np

from .._linalg import dot
from ..illuminants import whitepoints_cie1931
from ._color_space import ColorSpace


class PROLAB(ColorSpace):
    def __init__(self):
        super().__init__("proLab", ("L", "a", "b"), 0)
        # The matrix P in the article is
        #
        #  Q = ( Q 0 )
        #      ( q 1 )
        #
        self.Q = np.array(
            [
                [79.4725, 486.6610, 153.7311],
                [649.9038, -595.4477, -20.4498],
                [50.8625, 194.9377, -223.4334],
            ]
        )
        self.q = np.array([0.7947, 3.8666, 1.5373])
        self.Qinv = np.linalg.inv(self.Q)

    def from_xyz100(self, xyz):
        xyz = (xyz.T / whitepoints_cie1931["D65"]).T
        return dot(self.Q, xyz) / (dot(self.q, xyz) + 1)

    def to_xyz100(self, lab):
        y = dot(self.Qinv, lab)
        xyz = y / (1 - dot(self.q, y))
        xyz = (xyz.T * whitepoints_cie1931["D65"]).T
        return xyz
