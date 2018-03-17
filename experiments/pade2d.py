# -*- coding: utf-8 -*-
#
from __future__ import print_function, division

import numpy


def _print_triangle(alpha, degree):
    for d in range(degree+1):
        n = d*(d+1)//2
        print(alpha[n:n+d+1])
    return


def _create_triangle(alpha, degree):
    return [
        [alpha[d*(d+1)//2 + i] for i in range(d+1)]
        for d in range(degree+1)
        ]


def _evaluate_2d_polynomial(xy, alpha):
    x, y = xy
    out = 0.0
    for a in alpha:
        d = len(a)
        for i in range(d):
            out += a[i] * x**(d-i) * y**i
    return out



class Pade2d(object):
    '''Pad'e polynomial in from R^2 to R^2, i.e., both components are Pad'e
    functions in x and y. The function is characterized by four sets of
    coefficients, two for the components and two for numerator and denominator
    each.
    '''
    def __init__(self, degrees, shift):
        self.degrees = degrees
        self.shift = shift

        # Subtract 1 for each polynomial since the constant coefficient is fixed.
        self.num_coefficients = [(d+1)*(d+2)//2 - 1 for d in degrees]

        self.total_num_coefficients = sum(self.num_coefficients)

        # Choose the coefficiens to create the identity function
        self.alpha = numpy.zeros(self.total_num_coefficients)
        i0 = 0
        self.alpha[i0] = 1.0
        j0 = self.num_coefficients[0] + self.num_coefficients[1] + 1
        self.alpha[j0] = 1.0
        return

    def eval(self, xy):
        n = 0

        xy = (xy.T - self.shift).T

        alpha1 = numpy.concatenate([
            [0.0], self.alpha[n:n+self.num_coefficients[0]]
            ])
        a1 = _create_triangle(alpha1, self.degrees[0])
        n += self.num_coefficients[0]

        alpha2 = numpy.concatenate([
            [1.0], self.alpha[n:n+self.num_coefficients[1]]
            ])
        a2 = _create_triangle(alpha2, self.degrees[1])
        n += self.num_coefficients[1]

        beta1 = numpy.concatenate([
            [0.0], self.alpha[n:n+self.num_coefficients[2]]
            ])
        b1 = _create_triangle(beta1, self.degrees[2])
        n += self.num_coefficients[2]

        beta2 = numpy.concatenate([
            [1.0], self.alpha[n:n+self.num_coefficients[3]]
            ])
        b2 = _create_triangle(beta2, self.degrees[3])

        return numpy.array([
            _evaluate_2d_polynomial(xy, a1) / _evaluate_2d_polynomial(xy, a2),
            _evaluate_2d_polynomial(xy, b1) / _evaluate_2d_polynomial(xy, b2),
            ])

    def print(self):
        n = 0

        alpha1 = numpy.concatenate([
            [0.0], self.alpha[n:n+self.num_coefficients[0]]
            ])
        _print_triangle(alpha1, self.degrees[0])
        print()
        n += self.num_coefficients[0]

        alpha2 = numpy.concatenate([
            [1.0], self.alpha[n:n+self.num_coefficients[1]]
            ])
        _print_triangle(alpha2, self.degrees[1])
        print()
        n += self.num_coefficients[1]

        beta1 = numpy.concatenate([
            [0.0], self.alpha[n:n+self.num_coefficients[2]]
            ])
        _print_triangle(beta1, self.degrees[2])
        print()
        n += self.num_coefficients[2]

        beta2 = numpy.concatenate([
            [1.0], self.alpha[n:n+self.num_coefficients[3]]
            ])
        _print_triangle(beta2, self.degrees[3])
        print()
        return
