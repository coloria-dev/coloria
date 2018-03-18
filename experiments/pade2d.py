# -*- coding: utf-8 -*-
#
from __future__ import print_function, division

import numpy
import sympy


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
    def __init__(self, degrees, shift, alpha=None):
        self.degrees = degrees
        self.shift = shift

        # Subtract 1 for each polynomial since the constant coefficient is
        # fixed.
        self.num_coefficients = [(d+1)*(d+2)//2 - 1 for d in degrees]

        self.total_num_coefficients = sum(self.num_coefficients)

        if alpha is None:
            # Choose the coefficiens to create the identity function
            alpha = numpy.zeros(self.total_num_coefficients)
            i0 = 0
            alpha[i0] = 1.0
            j0 = self.num_coefficients[0] + self.num_coefficients[1] + 1
            alpha[j0] = 1.0

        self.set_alpha(alpha)
        return

    def set_alpha(self, alpha):
        assert len(alpha) == self.total_num_coefficients

        self.alpha = alpha

        n = 0

        alpha1 = numpy.concatenate([
            [0.0], alpha[n:n+self.num_coefficients[0]]
            ])
        self.a1 = _create_triangle(alpha1, self.degrees[0])
        n += self.num_coefficients[0]

        alpha2 = numpy.concatenate([
            [1.0], alpha[n:n+self.num_coefficients[1]]
            ])
        self.a2 = _create_triangle(alpha2, self.degrees[1])
        n += self.num_coefficients[1]

        beta1 = numpy.concatenate([
            [0.0], alpha[n:n+self.num_coefficients[2]]
            ])
        self.b1 = _create_triangle(beta1, self.degrees[2])
        n += self.num_coefficients[2]

        beta2 = numpy.concatenate([
            [1.0], alpha[n:n+self.num_coefficients[3]]
            ])
        self.b2 = _create_triangle(beta2, self.degrees[3])
        return

    def eval(self, xy):
        xy = (xy.T - self.shift).T

        p1 = _evaluate_2d_polynomial(xy, self.a1)
        q1 = _evaluate_2d_polynomial(xy, self.a2)

        p2 = _evaluate_2d_polynomial(xy, self.b1)
        q2 = _evaluate_2d_polynomial(xy, self.b2)

        return numpy.array([p1 / q1, p2 / q2])

    def jac(self, xy):
        '''Get the Jacobian at (x, y).
        '''
        x = sympy.Symbol('x')
        y = sympy.Symbol('y')

        return

    def print(self):
        for vals in [self.a1, self.a2, self.b1, self.b2]:
            for val in vals:
                print(val)
            print()
        return
