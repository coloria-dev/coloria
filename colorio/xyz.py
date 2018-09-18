# -*- coding: utf-8 -*-
#


# pylint: disable=no-self-use
class XYZ(object):
    def __init__(self):
        self.labels = ["X", "Y", "Z"]
        return

    def from_xyz100(self, xyz):
        return xyz / 100

    def to_xyz100(self, xyz):
        return xyz * 100
