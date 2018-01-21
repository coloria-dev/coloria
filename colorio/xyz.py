# -*- coding: utf-8 -*-
#
from __future__ import division

from .import srgb


def srgb_gamut(filename='srgb-xyz.vtu', n=50):
    srgb.show_gamut(filename, lambda x: x, n=n)
    return
