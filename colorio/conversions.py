# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from . import observers


def spectrum_to_xyz(spectrum, observer=observers.cie_1931_2):
    '''Computes the tristimulus values XYZ from a given spectrum for a given
    observer via

    X_i = int_lambda spectrum_i(lambda) * observer_i(lambda) dlambda.
    '''
    lambda_o, data_o = observer()

    # CIE 15:2004 says about the D65 illuminant:
    # If values at other wavelengths than printed in Table 1 of the standard
    # (CIE, 1998c) at 1 nm intervals are needed, linear interpolation should be
    # used.
    lambda_s, data_s = spectrum()

    # form the union of lambdas
    lmbda = numpy.sort(numpy.unique(numpy.concatenate([lambda_o, lambda_s])))

    # interpolate data
    idata_o = numpy.array([numpy.interp(lmbda, lambda_o, d) for d in data_o])
    idata_s = numpy.interp(lmbda, lambda_s, data_s)

    print(idata_o[0])

    exit(1)
    return


def xyz_to_xyy(xyz):
    return numpy.stack([xyz[:2] / numpy.sum(xyz, axis=0), xyz[1]])


def xyy_to_xyz(xyy):
    return numpy.stack([
        xyy[2] / xyy[1] * xyy[0],
        xyy[2],
        xyy[2] / xyy[1] * (1 - xyy[0] - xyy[1]),
        ])


def xyz_to_srgb1(xyz):
    # https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
    M = numpy.array([
        [+3.2406, -1.5372, -0.4986],
        [-0.9689, +1.8758, +0.0415],
        [+0.0557, -0.2040, +1.0570],
        ])
    srgb_linear = numpy.dot(M, xyz)
    a = 0.055
    srgb = numpy.array([
        12.92 * c if c <= 0.0031308 else
        (1+a) * c**(1/2.4) - a
        for c in srgb_linear
        ])
    return srgb


def srgb1_to_xyz(srgb1):
    # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
    a = 0.055
    srgb_linear = numpy.array([
        c / 12.92 if c <= 0.04045 else
        ((c+a) / (1+a)) ** 2.4
        for c in srgb1
        ])
    M = numpy.array([
        [0.4124, 0.3576, 0.1805],
        [0.2126, 0.7152, 0.0722],
        [0.0193, 0.1192, 0.9505],
        ])
    return numpy.dot(M, srgb_linear)


def srgb1_to_srgb256(srgb1):
    return (256 * srgb1).astype(int)


def srgb256_to_srgb1(srgb256):
    return srgb256 / 256
