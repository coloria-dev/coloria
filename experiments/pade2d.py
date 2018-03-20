# -*- coding: utf-8 -*-
#
from __future__ import print_function, division

import numpy
import sympy


def _create_tree(alpha, degree):
    return [
        numpy.array([alpha[d*(d+1)//2 + i] for i in range(d+1)])
        for d in range(degree+1)
        ]
    # return numpy.split(alpha, numpy.arange(1, degree+1))


def _get_xy_tree(xy, degree):
    '''Evaluates the entire tree of 2d mononomials.

    The return value is a list of arrays, where `out[k]` hosts the `2*k+1`
    values of the `k`th level of the tree

        (0, 0)
        (1, 0)   (0, 1)
        (2, 0)   (1, 1)   (0, 2)
          ...      ...      ...
    '''
    x, y = xy
    tree = [
        numpy.array([numpy.ones(x.shape, dtype=int)])
        ]
    for d in range(degree):
        tree.append(
            numpy.concatenate([tree[-1] * x, [tree[-1][-1]*y]])
            )
    return tree


def _get_dx_tree(xy, degree):
    '''
                                      0
                            1*(0, 0)  0
                  2*(1, 0)  1*(0, 1)  0
        3*(2, 0)  2*(1, 1)  1*(0, 2)  0
          ...       ...      ...     ...
    '''
    x, y = xy

    # build smaller tree
    one = numpy.array([numpy.ones(x.shape, dtype=int)])
    tree = [one]
    for d in range(1, degree):
        tree.append(
            numpy.concatenate([
                # Integer division `//` would be nice here, but
                # <https://github.com/sympy/sympy/issues/14542>.
                [tree[-1][0] / d * (d+1) * x],
                tree[-1] * y,
                ])
            )

    # append zeros
    zero = numpy.array([numpy.zeros(x.shape, dtype=int)])
    tree = [zero] + [numpy.concatenate([t, zero]) for t in tree]
    return tree


def _get_dy_tree(xy, degree):
    '''
        0
        0  1*(0, 0)
        0  1*(1, 0)  2*(0, 1)
        0  1*(2, 0)  2*(1, 1)  3*(0, 2)
       ...   ...       ...       ...
    '''
    x, y = xy

    one = numpy.array([numpy.ones(x.shape, dtype=int)])
    tree = [one]
    for d in range(1, degree):
        tree.append(
            numpy.concatenate([
                tree[-1] * x,
                # Integer division `//` would be nice here, but
                # <https://github.com/sympy/sympy/issues/14542>.
                [tree[-1][-1] / d * (d+1) * y]
                ])
            )

    # prepend zeros
    zero = numpy.array([numpy.zeros(x.shape, dtype=int)])
    tree = [zero] + [numpy.concatenate([zero, t]) for t in tree]
    return tree


def _eval_tree(xy_tree, coeff_tree):
    # print(xy_tree[1].shape, coeff_tree[1].shape)
    return numpy.sum([
        numpy.sum((xy_tree[k].T * coeff_tree[k]).T, axis=0)
        for k in range(len(coeff_tree))
        ], axis=0)


class Pade2d(object):
    '''Pad'e polynomial in from R^2 to R^2, i.e., both components are Pad'e
    functions in x and y. The function is characterized by four sets of
    coefficients, two for the components and two for numerator and denominator
    each.
    '''
    def __init__(self, degrees, alpha=None):
        self.degrees = degrees

        # Subtract 1 for each denominator polynomial since the constant
        # coefficient is fixed to 1.0.
        self.num_coefficients = [
            (degrees[0]+1)*(degrees[0]+2)//2,
            (degrees[1]+1)*(degrees[1]+2)//2 - 1,
            (degrees[2]+1)*(degrees[2]+2)//2,
            (degrees[3]+1)*(degrees[3]+2)//2 - 1,
            ]

        self.num_coefficients = [(d+1)*(d+2)//2 for d in degrees]

        # Subtract 1 for each denominator polynomial since the constant
        # coefficient is fixed to 1.0.
        self.total_num_coefficients = sum(self.num_coefficients) - 2

        if alpha is None:
            # Choose the coefficiens to create the identity function
            coeffs_ax = numpy.zeros(self.num_coefficients[0])
            coeffs_ax[1] = 1
            coeffs_bx = numpy.zeros(self.num_coefficients[1] - 1)

            coeffs_ay = numpy.zeros(self.num_coefficients[2])
            coeffs_ay[2] = 1
            coeffs_by = numpy.zeros(self.num_coefficients[3] - 1)

            alpha = numpy.concatenate([
                coeffs_ax, coeffs_bx, coeffs_ay, coeffs_by
                ])

        self.set_alpha(alpha)
        return

    def set_xy(self, xy):
        self.xy = xy
        self.xy_tree = _get_xy_tree(xy, max(self.degrees))
        self.dx_tree = _get_dx_tree(xy, max(self.degrees))
        self.dy_tree = _get_dy_tree(xy, max(self.degrees))
        return

    def set_alpha(self, alpha):
        assert len(alpha) == self.total_num_coefficients

        self.alpha = alpha

        num_coefficients = [(d+1)*(d+2)//2 for d in self.degrees]
        num_coefficients[1] -= 1
        num_coefficients[3] -= 1

        ax, bx, ay, by = numpy.split(alpha, numpy.cumsum(num_coefficients[:-1]))
        bx = numpy.concatenate([[1.0], bx])
        by = numpy.concatenate([[1.0], by])

        self.tree_ax = _create_tree(ax, self.degrees[0])
        self.tree_ay = _create_tree(ay, self.degrees[1])
        self.tree_bx = _create_tree(bx, self.degrees[2])
        self.tree_by = _create_tree(by, self.degrees[3])
        return

    def eval(self, xy=None):
        if xy is not None:
            self.set_xy(xy)

        ux = _eval_tree(self.xy_tree, self.tree_ax)
        vx = _eval_tree(self.xy_tree, self.tree_bx)
        uy = _eval_tree(self.xy_tree, self.tree_ay)
        vy = _eval_tree(self.xy_tree, self.tree_by)

        return numpy.array([ux / vx, uy / vy])

    def jac(self, xy=None):
        '''Get the Jacobian at (x, y).
        '''
        if xy is not None:
            self.set_xy(xy)

        ux = _eval_tree(self.xy_tree, self.tree_ax)
        vx = _eval_tree(self.xy_tree, self.tree_bx)
        uy = _eval_tree(self.xy_tree, self.tree_ay)
        vy = _eval_tree(self.xy_tree, self.tree_by)

        ux_dx = _eval_tree(self.dx_tree, self.tree_ax)
        vx_dx = _eval_tree(self.dx_tree, self.tree_bx)
        uy_dx = _eval_tree(self.dx_tree, self.tree_ay)
        vy_dx = _eval_tree(self.dx_tree, self.tree_by)

        ux_dy = _eval_tree(self.dy_tree, self.tree_ax)
        vx_dy = _eval_tree(self.dy_tree, self.tree_bx)
        uy_dy = _eval_tree(self.dy_tree, self.tree_ay)
        vy_dy = _eval_tree(self.dy_tree, self.tree_by)

        jac = numpy.array([
            [
                (ux_dx * vx - vx_dx * ux) / vx**2,
                (ux_dy * vx - vx_dy * ux) / vx**2,
            ],
            [
                (uy_dx * vy - vy_dx * uy) / vy**2,
                (uy_dy * vy - vy_dy * uy) / vy**2,
            ],
            ])
        return jac

    def print(self):
        for k, vals in enumerate([self.tree_ax, self.tree_bx, self.tree_ay, self.tree_by]):
            for val in vals:
                print(val)
            print()
        return
