from ._srgb import SrgbLinear


class ColorSpace:
    is_origin_well_defined = True
    srgb_linear = SrgbLinear()

    def __repr__(self):
        return f"<colorio color space {self.name}>"

    def to_xyz100(self, _):
        raise NotImplementedError("ColorSpace needs to implement to_xyz100()")

    def from_xyz100(self, _):
        raise NotImplementedError("ColorSpace needs to implement from_xyz100()")

    def to_rgb_linear(self, cs_coords):
        return self.srgb_linear.from_xyz100(self.to_xyz100(cs_coords))

    def to_rgb1(self, cs_coords):
        return self.srgb_linear.to_rgb1(self.to_rgb_linear(cs_coords))

    def to_rgb255(self, cs_coords):
        return self.srgb_linear.to_rgb255(self.to_rgb_linear(cs_coords))

    def from_rgb_linear(self, rgb_lin):
        return self.from_xyz100(self.srgb_linear.to_xyz100(rgb_lin))

    def from_rgb1(self, rgb1):
        return self.from_rgb_linear(self.srgb_linear.from_rgb1(rgb1))

    def from_rgb255(self, rgb255):
        return self.from_rgb_linear(self.srgb_linear.from_rgb255(rgb255))
