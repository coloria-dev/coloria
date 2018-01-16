# -*- coding: utf-8 -*-
#
import colorio


def test_dsad():
    xyz = [0.5, 0.4, 0.3]
    srgb = colorio.xyz_to_srgb(xyz)
    print(srgb)
    return
