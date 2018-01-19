# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy


def from_xyz(xyz):
    return numpy.concatenate([xyz[:2] / numpy.sum(xyz, axis=0), [xyz[1]]])


def to_xyz(xyy):
    return numpy.stack([
        xyy[2] / xyy[1] * xyy[0],
        xyy[2],
        xyy[2] / xyy[1] * (1 - xyy[0] - xyy[1]),
        ])
