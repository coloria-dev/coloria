import numpy

from ._color_space import ColorSpace


class XYY(ColorSpace):
    def __init__(self):
        super().__init__()
        self.labels = ["x", "y", "Y"]
        self.k0 = 2  # which index corresponds to luminosity

    def from_xyz100(self, xyz100):
        xyz = xyz100 / 100
        sum_xyz = numpy.sum(xyz, axis=0)
        x = xyz[0]
        y = xyz[1]
        return numpy.array([x / sum_xyz, y / sum_xyz, y])

    def to_xyz100(self, xyy):
        x, y, Y = xyy
        return numpy.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100
