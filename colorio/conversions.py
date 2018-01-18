# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from . import observers


def spectrum_to_xyz(spectrum, observer=observers.cie_1931_2()):
    '''Computes the tristimulus values XYZ from a given spectrum for a given
    observer via

    X_i = int_lambda spectrum_i(lambda) * observer_i(lambda) dlambda.

    In section 7, the technical report CIE Standard Illuminants for
    Colorimetry, 1999, gives a recommendation on how to perform the
    computation.
    '''
    lambda_o, data_o = observer
    lambda_s, data_s = spectrum

    # form the union of lambdas
    lmbda = numpy.sort(numpy.unique(numpy.concatenate([lambda_o, lambda_s])))

    # The technical document prescribes that the integration be performed "over
    # the wavelength range corresponding to the entire visible spectrum, 360 nm
    # to 830 nm.
    assert lmbda[0] < 361e-9
    assert lmbda[-1] > 829e-9

    # interpolate data
    idata_o = numpy.array([numpy.interp(lmbda, lambda_o, d) for d in data_o])
    # The technical report specifies the interpolation techniques, too:
    # ```
    # Use one of the four following methods to calculate needed but unmeasured
    # values of phi(l), R(l) or tau(l) within the range of measurements:
    #   1) the third-order polynomial interpolation (Lagrange) from the four
    #      neighbouring data points around the point to be interpolated, or
    #   2) cubic spline interpolation formula, or
    #   3) a fifth order polynomial interpolation formula from the six
    #      neighboring data points around the point to be interpolated, or
    #   4) a Sprague interpolation (see Seve, 2003).
    # ```
    # Well, don't do that but simply use linear interpolation now. We only use
    # the midpoint rule for integration anyways.
    idata_s = numpy.interp(lmbda, lambda_s, data_s)

    values = numpy.dot(idata_o, idata_s)

    # scale the values such that Y=1
    values /= values[1]

    return values


def xyz_to_xyy(xyz):
    return numpy.concatenate([xyz[:2] / numpy.sum(xyz, axis=0), [xyz[1]]])


def xyy_to_xyz(xyy):
    return numpy.stack([
        xyy[2] / xyy[1] * xyy[0],
        xyy[2],
        xyy[2] / xyy[1] * (1 - xyy[0] - xyy[1]),
        ])


def xyz_to_srgb1(xyz):
    # https://en.wikipedia.org/wiki/SRGB#The_forward_transformation_(CIE_XYZ_to_sRGB)
    # http://www.color.org/srgb.pdf
    M = numpy.array([
        [+3.2406255, -1.537208, -0.4986286],
        [-0.9689307, +1.8757561, +0.0415175],
        [+0.0557101, -0.2040211, +1.0569959],
        ])
    srgb_linear = numpy.dot(M, xyz)

    a = 0.055
    is_smaller = srgb_linear <= 0.0031308
    is_greater = numpy.logical_not(is_smaller)

    srgb = srgb_linear
    srgb[is_smaller] *= 12.92
    srgb[is_greater] = (1+a) * srgb[is_greater]**(1/2.4) - a
    return srgb


def srgb1_to_xyz(srgb1):
    srgb_linear = numpy.array(srgb1, dtype=float)

    # https://en.wikipedia.org/wiki/SRGB#The_reverse_transformation
    is_smaller = srgb_linear <= 0.04045
    is_greater = numpy.logical_not(is_smaller)

    a = 0.055
    srgb_linear[is_smaller] /= 12.92
    srgb_linear[is_greater] = ((srgb_linear[is_greater] + a) / (1+a))**2.4

    M = numpy.array([
        [+3.2406255, -1.537208, -0.4986286],
        [-0.9689307, +1.8757561, +0.0415175],
        [+0.0557101, -0.2040211, +1.0569959],
        ])
    return numpy.linalg.solve(M, srgb_linear)


def srgb1_to_srgb256(srgb1):
    return (256 * srgb1).astype(int)


def srgb256_to_srgb1(srgb256):
    return srgb256 / 256
