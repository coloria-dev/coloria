# https://bottosson.github.io/posts/oklab/
import numpy

from .._linalg import dot
from ._color_space import ColorSpace


class OKLAB(ColorSpace):
    def __init__(self):
        super().__init__("Oklab", ("L", "a", "b"), 0)
        self.M1 = numpy.array(
            [
                [0.8189330101, 0.3618667424, -0.1288597137],
                [0.0329845436, 0.9293118715, 0.0361456387],
                [0.0482003018, 0.2643662691, 0.6338517070],
            ]
        )
        self.M1inv = numpy.linalg.inv(self.M1)
        self.M2 = numpy.array(
            [
                [0.2104542553, +0.7936177850, -0.0040720468],
                [+1.9779984951, -2.4285922050, +0.4505937099],
                [+0.0259040371, +0.7827717662, -0.8086757660],
            ]
        )
        self.M2inv = numpy.linalg.inv(self.M2)

    def from_xyz100(self, xyz):
        return dot(self.M2, numpy.cbrt(dot(self.M1, xyz)))

    def to_xyz100(self, xyz):
        return dot(self.M1inv, dot(self.M2inv, xyz) ** 3)
