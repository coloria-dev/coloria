"""
https://bottosson.github.io/posts/oklab/
"""
import npx
import numpy as np
from numpy.typing import ArrayLike

from ._color_space import ColorSpace
from ._helpers import register


class OKLAB(ColorSpace):
    name = "Oklab"
    labels = ("L", "a", "b")
    k0 = 0

    def __init__(self):
        self.M1 = np.array(
            [
                [0.8189330101, 0.3618667424, -0.1288597137],
                [0.0329845436, 0.9293118715, 0.0361456387],
                [0.0482003018, 0.2643662691, 0.6338517070],
            ]
        )
        self.M1inv = np.linalg.inv(self.M1)
        self.M2 = np.array(
            [
                [0.2104542553, +0.7936177850, -0.0040720468],
                [+1.9779984951, -2.4285922050, +0.4505937099],
                [+0.0259040371, +0.7827717662, -0.8086757660],
            ]
        )
        self.M2inv = np.linalg.inv(self.M2)

    def from_xyz100(self, xyz100: ArrayLike) -> np.ndarray:
        xyz = np.asarray(xyz100) / 100
        return npx.dot(self.M2, np.cbrt(npx.dot(self.M1, xyz)))

    def to_xyz100(self, lab: ArrayLike) -> np.ndarray:
        return npx.dot(self.M1inv, npx.dot(self.M2inv, lab) ** 3) * 100


register("oklab", OKLAB())
