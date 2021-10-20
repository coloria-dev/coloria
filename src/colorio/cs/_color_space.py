class ColorSpace:
    name = "unknown"
    labels = ("", "", "")
    is_origin_well_defined = True
    k0 = None

    def __repr__(self):
        return f"<colorio color space {self.name}>"

    def to_xyz100(self, _):
        raise NotImplementedError("ColorSpace needs to implement to_xyz100()")

    def from_xyz100(self, _):
        raise NotImplementedError("ColorSpace needs to implement from_xyz100()")

    @property
    def lightness_label(self):
        assert self.k0 is not None
        return self.labels[self.k0]

    @property
    def hue_labels(self):
        assert self.k0 is not None
        hue_idx = [True, True, True]
        hue_idx[self.k0] = False
        return [label for k, label in enumerate(self.labels) if k != self.k0]
