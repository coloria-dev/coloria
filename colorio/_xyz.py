from ._color_space import ColorSpace


class XYZ(ColorSpace):
    def __init__(self):
        self.labels = ["X", "Y", "Z"]
        return

    def from_xyz100(self, xyz):
        return xyz / 100

    def to_xyz100(self, xyz):
        return xyz * 100
