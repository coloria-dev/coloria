# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .linalg import dot, solve


class ICtCp(object):
    '''
    ICtCp color model.
    <https://en.wikipedia.org/wiki/ICtCp>
    '''
    def __init__(self):
        self.M1 = numpy.array([
            [1688, 2146, 262],
            [683, 2951, 462],
            [99, 309, 3688]
            ]) / 4096

        # From <https://doi.org/10.5594/SMPTE.ST2084.2014>
        self.m1 = 2610 / 4096 / 4
        self.m2 = 2523 / 4096 / 128
        self.c1 = 3424 / 4096
        self.c2 = 2413 / 4096 * 32
        self.c3 = 2392 / 4096 * 32

        self.M2 = numpy.array([
            [2048, 2048, 0],
            [6610, -13613, 7003],
            [17933, -17390, -543],
            ]) / 4096
        return

    def from_rec2100(self, rgb):
        lms = dot(self.M1, rgb)

        lms_ = (
            (self.c1 + self.c2 * lms**self.m1) / (1 + self.c3 * lms**self.m1)
            )**self.m2

        ictcp = dot(self.M2, lms_)
        return ictcp

    def to_rec2100(self, ictcp):
        lms_ = solve(self.M2, ictcp)

        t = lms_**(1/self.m2) - self.c1
        # This next line is part of the model, but really it shouldn't occur
        # for sane input data.
        # t[t < 0] = 0.0
        lms = (t / (self.c2 - self.c3*lms_**(1/self.m2)))**(1/self.m1)

        rgb = solve(self.M1, lms)
        return rgb
