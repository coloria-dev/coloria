import numpy

from .._exceptions import ColorioError
from ._color_space import ColorSpace


class XYY(ColorSpace):
    """xyY colorspace with Y scaled from 0 to 1."""

    def __init__(self, Y_scaling: int):
        if not Y_scaling in [1, 100]:
            raise ColorioError("Y_scaling needs to be 1 or 100.")

        super().__init__(
            f"xyY_{Y_scaling}",
            ("x", "y", f"Y{Y_scaling}"),
            2,
            is_origin_well_defined=False,
        )
        self.Y_scaling = Y_scaling

    def from_xyz100(self, xyz100):
        if numpy.any(xyz100 < 0):
            raise ColorioError("Negative XYZ100 value.")

        xyz = numpy.asarray(xyz100)
        if self.Y_scaling == 1:
            xyz = xyz100 / 100
        sum_xyz = numpy.sum(xyz, axis=0)
        x = xyz[0]
        y = xyz[1]
        return numpy.array([x / sum_xyz, y / sum_xyz, y])

    def to_xyz100(self, xyy):
        if numpy.any(xyy < 0):
            raise ColorioError("Negative xyY value.")
        x, y, Y = xyy
        out = numpy.array([Y / y * x, Y, Y / y * (1 - x - y)])
        if self.Y_scaling == 1:
            out *= 100
        return out
