# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from . import conversions
from . import observers


def white_point(illuminant, observer=observers.cie_1931_2()):
    '''From <https://en.wikipedia.org/wiki/White_point>:
    The white point of an illuminant is the chromaticity of a white object
    under the illuminant.
    '''
    return conversions.spectrum_to_xyz(illuminant, observer)


def a(interval=1):
    '''CIE Standard Illuminants for Colorimetry, 1999:
    CIE standard illuminant A is intended to represent typical, domestic,
    tungsten-filament lighting. Its relative spectral power distribution is
    that of a Planckian radiator at a temperature of approximately 2856 K. CIE
    standard illuminant A should be used in all applications of colorimetry
    involving the use of incandescent lighting, unless there are specific
    reasons for using a different illuminant.
    '''
    # https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_A
    lmbda = numpy.arange(300, 831, interval)
    c2 = 1.435e7
    color_temp = 2848
    vals = 100 * (560 / lmbda)**5 * (
        (numpy.exp(c2 / (color_temp * 560)) - 1) /
        (numpy.exp(c2 / (color_temp * lmbda)) - 1)
        )
    return lmbda, vals


def d(nominal_temperature):
    '''CIE D-series illuminants.

    The technical report `Colorimetry, 3rd edition, 2004` gives the data for
    D50, D55, and D65 explicitly, but also explains how it's computes for S0,
    S1, S2. Values are given at 5nm resolution in the document, but really
    every other value is just interpolated. Hence, only provide 10 nm data
    here.
    '''
    # From CIE 15:2004. Colorimetry, 3rd edition, 2004 (page 69, note 5):
    #
    # The method required to calculate the values for the relative spectral
    # power distributions of illuminants D50, D55, D65, and D75, in Table T.1
    # is as follows
    #   1. Multiply the nominal correlated colour temperature (5000 K, 5500 K,
    #      6500 K or 7500 K) by 1,4388/1,4380.
    #   2. Calculate XD and YD using the equations given in the text.
    #   3. Calculate M1 and M2 using the equations given in the text.
    #   4. Round M1 and M2 to three decimal places.
    #   5. Calculate S(lambda) every 10 nm by
    #        S(lambda) = S0(lambda) + M1 S1(lambda) + M2 S2(lambda)
    #      using values of S0(lambda), S1(lambda) and S2(lambda) from
    #      Table T.2.
    #   6. Interpolate the 10 nm values of S(lambda) linearly to obtain values
    #      at intermediate wavelengths.
    tcp = 1.4388/1.4380 * nominal_temperature

    if 4000 <= tcp <= 7000:
        xd = -4.6070e9/tcp**3 + 2.9678e6/tcp**2 + 0.09911e3/tcp + 0.244063
    else:
        assert 7000 < tcp <= 25000
        xd = -2.0064e9/tcp**3 + 1.9018e6/tcp**2 + 0.24748e3/tcp + 0.237040

    yd = -3.000*xd**2 + 2.870*xd - 0.275

    m1 = (-1.3515 - 1.7703*xd + 5.9114*yd) / (0.0241 + 0.2562*xd - 0.7341*yd)
    m2 = (+0.0300 - 31.4424*xd + 30.0717*yd) / (0.0241 + 0.2562*xd - 0.7341*yd)

    m1 = numpy.around(m1, decimals=3)
    m2 = numpy.around(m2, decimals=3)

    lmbda = numpy.arange(300, 831, 10)
    # The standard gives values every 5nm, but every other value is just an
    # interpolation from 10nm data. No idea why they are listed in the standard
    s = numpy.array([
        # 300 nm:
        [0.04, 0.02, 0.00],
        [6.00, 4.50, 2.00],
        [29.60, 22.40, 4.00],
        [55.30, 42.00, 8.50],
        [57.30, 40.60, 7.80],
        [61.80, 41.60, 6.70],
        [61.50, 38.00, 5.30],
        [68.80, 42.40, 6.10],
        [63.40, 38.50, 3.00],
        [65.80, 35.00, 1.20],
        # 400 nm:
        [94.80, 43.40, -1.10],
        [104.80, 46.30, -0.50],
        [105.90, 43.90, -0.70],
        [96.80, 37.10, -1.20],
        [113.90, 36.70, -2.60],
        [125.60, 35.90, -2.90],
        [125.50, 32.60, -2.80],
        [121.30, 27.90, -2.60],
        [121.30, 24.30, -2.60],
        [113.50, 20.10, -1.80],
        # 500 nm:
        [113.10, 16.20, -1.50],
        [110.80, 13.20, -1.30],
        [106.50, 8.60, -1.20],
        [108.80, 6.10, -1.00],
        [105.30, 4.20, -0.50],
        [104.40, 1.90, -0.30],
        [100.00, 0.00, 0.00],
        [96.00, -1.60, 0.20],
        [95.10, -3.50, 0.50],
        [89.10, -3.50, 2.10],
        # 600 nm:
        [90.50, -5.80, 3.20],
        [90.30, -7.20, 4.10],
        [88.40, -8.60, 4.70],
        [84.00, -9.50, 5.10],
        [85.10, -10.90, 6.70],
        [81.90, -10.70, 7.30],
        [82.60, -12.00, 8.60],
        [84.90, -14.00, 9.80],
        [81.30, -13.60, 10.20],
        [71.90, -12.00, 8.30],
        # 700 nm:
        [74.30, -13.30, 9.60],
        [76.40, -12.90, 8.50],
        [63.30, -10.60, 7.00],
        [71.70, -11.60, 7.60],
        [77.00, -12.20, 8.00],
        [65.20, -10.20, 6.70],
        [47.70, -7.80, 5.20],
        [68.60, -11.20, 7.40],
        [65.00, -10.40, 6.80],
        [66.00, -10.60, 7.00],
        # 800 nm:
        [61.00, -9.70, 6.40],
        [53.30, -8.30, 5.50],
        [58.90, -9.30, 6.10],
        [61.90, -9.80, 6.50],
        ]).T

    return lmbda, s[0] + m1*s[1] + m2*s[2]


def d50():
    '''CIE illuminant D50, mid-morning/mid-afternoon daylight, at 10nm
    resolution.
    '''
    return d(5000)


def d55():
    '''CIE illuminant D55, mid-morning/mid-afternoon daylight, at 10nm
    resolution.
    '''
    return d(5500)


def d65():
    '''CIE standard illuminant D65, sampled at 10nm intervals.
    '''
    return d(6500)


def d75():
    '''CIE illuminant D75
    '''
    return d(7500)


def e():
    '''This is a hypothetical reference radiator. All wavelengths in CIE
    illuminant E are weighted equally with a relative spectral power of 100.0.
    '''
    lmbda = numpy.arange(300, 831)
    data = numpy.full(lmbda.shape, 100.0)
    return lmbda, data


def f2():
    # http://www.npsg.uwaterloo.ca/data/illuminant.php
    lmbda = numpy.arange(300, 781, 5)
    vals = numpy.array([
        1.18, 1.48, 1.84, 2.15, 3.44, 15.69, 3.85, 3.74, 4.19, 4.62, 5.06,
        34.98, 11.81, 6.27, 6.63, 6.93, 7.19, 7.40, 7.54, 7.62, 7.65, 7.62,
        7.62, 7.45, 7.28, 7.15, 7.05, 7.04, 7.16, 7.47, 8.04, 8.88, 10.01,
        24.88, 16.64, 14.59, 16.16, 17.56, 18.62, 21.47, 22.79, 19.29, 18.66,
        17.73, 16.54, 15.21, 13.80, 12.36, 10.95, 9.65, 8.40, 7.32, 6.31, 5.43,
        4.68, 4.02, 3.45, 2.96, 2.55, 2.19, 1.89, 1.64, 1.53, 1.27, 1.10, 0.99,
        0.88, 0.76, 0.68, 0.61, 0.56, 0.54, 0.51, 0.47, 0.47, 0.43, 0.46, 0.47,
        0.4, 0.33, 0.27,
        ])
    return lmbda, vals
