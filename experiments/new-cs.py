# -*- coding: utf-8 -*-
#
from __future__ import print_function, division

import os

import matplotlib.pyplot as plt
import numpy
from scipy.optimize import leastsq, least_squares
import yaml

import colorio


def evaluate_2d_polynomial(xy, alpha):
    x, y = xy
    out = 0.0
    for a in alpha:
        d = len(a)
        for i in range(d):
            out += a[i] * x**(d-i) * y**i
    return out


def create_triangle(alpha, degree):
    return [
        [alpha[d*(d+1)//2 + i] for i in range(d+1)]
        for d in range(degree+1)
        ]


def transform(xy, shift, alpha, num_coefficients, poly_degrees):
    n = 0

    xy = (xy.T - shift).T

    alpha1 = numpy.concatenate([[0.0], alpha[n:n+num_coefficients[0]]])
    a1 = create_triangle(alpha1, poly_degrees[0])
    n += num_coefficients[0]

    alpha2 = numpy.concatenate([[1.0], alpha[n:n+num_coefficients[1]]])
    a2 = create_triangle(alpha2, poly_degrees[1])
    n += num_coefficients[1]

    beta1 = numpy.concatenate([[0.0], alpha[n:n+num_coefficients[2]]])
    b1 = create_triangle(beta1, poly_degrees[2])
    n += num_coefficients[2]

    beta2 = numpy.concatenate([[1.0], alpha[n:n+num_coefficients[3]]])
    b2 = create_triangle(beta2, poly_degrees[3])
    n += num_coefficients[3]

    return numpy.array([
        evaluate_2d_polynomial(xy, a1) / evaluate_2d_polynomial(xy, a2),
        evaluate_2d_polynomial(xy, b1) / evaluate_2d_polynomial(xy, b2),
        ])


# pylint: disable=too-many-arguments
def get_ecc(alpha, num_coefficients, poly_degrees, centers, points, shift):
    '''Get eccentricities of ellipses.
    '''
    A = []
    B = []
    for center, pts in zip(centers, points):
        tcenter = transform(center, shift, alpha, num_coefficients, poly_degrees)
        X = (
            transform(pts, shift, alpha, num_coefficients, poly_degrees).T
            - tcenter
            ).T

        def f_ellipse(a_b_theta, x=X):
            a, b, theta = a_b_theta
            return (
                + a**2 * (x[0] * numpy.cos(theta) + x[1] * numpy.sin(theta))**2
                + b**2 * (x[0] * numpy.sin(theta) - x[1] * numpy.cos(theta))**2
                - 1.0
                )

        def jac_ellipse(a_b_theta, x=X):
            a, b, theta = a_b_theta
            return numpy.array([
                + 2*a * (x[0] * numpy.cos(theta) + x[1] * numpy.sin(theta))**2,
                + 2*b * (x[0] * numpy.sin(theta) - x[1] * numpy.cos(theta))**2,
                + a**2 * 2*(x[0] * numpy.cos(theta) + x[1] * numpy.sin(theta))
                * (-x[0] * numpy.sin(theta) + x[1] * numpy.cos(theta))
                + b**2 * 2*(x[0] * numpy.sin(theta) - x[1] * numpy.cos(theta))
                * (x[0] * numpy.cos(theta) + x[1] * numpy.sin(theta)),
                ]).T

        (a, b, _), _ = \
            leastsq(f_ellipse, [1.0, 1.0, 0.0], Dfun=jac_ellipse)

        A.append(1/a)
        B.append(1/b)

        # from matplotlib.patches import Ellipse
        # plt.plot(*X, 'x')
        # ax = plt.gca()
        # e = Ellipse(
        #     xy=[0.0, 0.0],
        #     width=2/a, height=2/b,
        #     angle=theta / numpy.pi * 180
        #     )
        # ax.add_artist(e)
        # e.set_alpha(0.5)
        # e.set_facecolor('k')
        # plt.show()

    return numpy.concatenate([A, B])
    # return numpy.log(ab / target_radius)


def print_triangle(alpha, degree):
    for d in range(degree+1):
        n = d*(d+1)//2
        print(alpha[n:n+d+1])
    return


def print_parameters(alpha, num_coefficients, poly_degrees):
    n = 0

    alpha1 = numpy.concatenate([[0.0], alpha[n:n+num_coefficients[0]]])
    print_triangle(alpha1, poly_degrees[0])
    print()
    n += num_coefficients[0]

    alpha2 = numpy.concatenate([[1.0], alpha[n:n+num_coefficients[1]]])
    print_triangle(alpha2, poly_degrees[1])
    print()
    n += num_coefficients[1]

    beta1 = numpy.concatenate([[0.0], alpha[n:n+num_coefficients[2]]])
    print_triangle(beta1, poly_degrees[2])
    print()
    n += num_coefficients[2]

    beta2 = numpy.concatenate([[1.0], alpha[n:n+num_coefficients[3]]])
    print_triangle(beta2, poly_degrees[3])
    print()
    n += num_coefficients[3]
    return


def _main():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, '../colorio/data/macadam1942/table3.yaml')) as f:
        data = yaml.safe_load(f)

    centers = []
    points = []
    for datak in data:
        # collect ellipse points
        _, _, _, _, delta_y_delta_x, delta_s = numpy.array(datak['data']).T
        if len(delta_s) < 2:
            continue
        centers.append([datak['x'], datak['y']])
        points.append((
            (
                numpy.array([numpy.ones(delta_y_delta_x.shape[0]), delta_y_delta_x])
                / numpy.sqrt(1 + delta_y_delta_x**2) * delta_s
            ).T + centers[-1]
            ).T)


    centers = numpy.array(centers)

    # shift white point to center (0.0, 0.0)
    whitepoint = [1/3, 1/3]

    poly_degrees = [6, 6, 6, 6]

    # Subtract 1 for each polynomial since the constant coefficient is fixed.
    num_coefficients = [(d+1)*(d+2)//2 - 1 for d in poly_degrees]

    def f2(alpha):
        ecc = get_ecc(
            alpha, num_coefficients, poly_degrees, centers, points, whitepoint
            )
        # compute standard deviation of ecc
        average = numpy.sum(ecc) / len(ecc)
        return (ecc - average) / average

    # Create the identity function as initial guess
    total_num_parameters = numpy.sum(num_coefficients)
    print('num parameters: {}'.format(total_num_parameters))
    coeff0 = numpy.zeros(total_num_parameters)
    i0 = 0
    coeff0[i0] = 1.0
    j0 = num_coefficients[0] + num_coefficients[1] + 1
    coeff0[j0] = 1.0
    print('\ninitial parameters:')
    print_parameters(coeff0, num_coefficients, poly_degrees)
    # out = leastsq(f, coeff0, full_output=True)
    # print(out)
    # exit(1)
    # coeff1, _ = leastsq(f2, coeff0, maxfev=10000)

    # Levenberg-Marquardt (lm) is better suited for small, dense, unconstrained
    # problems, but it needs more conditions than parameters.
    out = least_squares(f2, coeff0, method='trf')
    coeff1 = out.x
    print('\noptimal parameters:')
    print_parameters(coeff1, num_coefficients, poly_degrees)

    # plot statistics
    ecc0 = get_ecc(coeff0, num_coefficients, poly_degrees, centers, points, whitepoint)
    ecc1 = get_ecc(coeff1, num_coefficients, poly_degrees, centers, points, whitepoint)
    plt.plot(ecc0, label='ecc before')
    plt.plot(ecc1, label='ecc opt')
    plt.legend()

    # Plot unperturbed MacAdam
    plt.figure()
    colorio.plot_macadam(
        scaling=10,
        plot_standard_deviations=True
        )

    # Plot perturbed MacAdam
    plt.figure()
    colorio.plot_macadam(
        scaling=10,
        xy_to_2d=lambda xy: transform(xy, whitepoint, coeff1, num_coefficients, poly_degrees),
        plot_standard_deviations=True
        )

    plt.show()
    return


if __name__ == '__main__':
    _main()
