# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy


# pylint: disable=no-self-use
class XYY(object):
    def __init__(self):
        return

    def from_xyz100(self, xyz):
        sum_xyz = numpy.sum(xyz, axis=0)
        x, y, _ = xyz
        return numpy.array([x/sum_xyz, y/sum_xyz, y])

    def to_xyz100(self, xyy):
        x, y, Y = xyy
        return numpy.array([Y/y*x, Y, Y/y * (1-x-y)])
