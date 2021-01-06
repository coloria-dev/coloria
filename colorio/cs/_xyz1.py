from ._color_space import ColorSpace


class XYZ1(ColorSpace):
    """XYZ colorspace with the values scaled from 0 to 1."""

    def __init__(self):
        super().__init__("XYZ", ("X", "Y", "Z"), None)

    def from_xyz100(self, xyz):
        return xyz / 100

    def to_xyz100(self, xyz):
        return xyz * 100
