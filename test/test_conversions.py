# -*- coding: utf-8 -*-
#
import colorio


def test_spectrum_to_xyz():
    print(colorio.spectrum_to_xyz(colorio.illuminants.d65()))
    return


if __name__ == '__main__':
    test_spectrum_to_xyz()
