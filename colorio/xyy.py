# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy


class XYY(object):
    def __init__(self):
        self.labels = ["x", "y", "Y"]
        return

    def from_xyz100(self, xyz):
        sum_xyz = numpy.sum(xyz, axis=0)
        x = xyz[0]
        y = xyz[1]
        return numpy.array([x / sum_xyz, y / sum_xyz, y / 100])

    def to_xyz100(self, xyy):
        x, y, Y = xyy
        return numpy.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100
