# -*- coding: utf-8 -*-
#
from __future__ import print_function, division

import os

import matplotlib.pyplot as plt
import numpy
from scipy.optimize import leastsq, least_squares
import yaml

import colorio

from pade2d import Pade2d


def f_ellipse(a_b_theta, x):
    a, b, theta = a_b_theta
    cos = numpy.cos(theta)
    sin = numpy.sin(theta)
    return (
        + a**2 * (x[0]*cos + x[1]*sin)**2
        + b**2 * (x[0]*sin - x[1]*cos)**2
        - 1.0
        )


def jac_ellipse(a_b_theta, x):
    a, b, theta = a_b_theta
    cos = numpy.cos(theta)
    sin = numpy.sin(theta)
    return numpy.array([
        + 2*a * (x[0]*cos + x[1]*sin)**2,
        #
        + 2*b * (x[0]*sin - x[1]*cos)**2,
        #
        + a**2 * 2*(x[0]*cos + x[1]*sin) * (-x[0]*sin + x[1]*cos)
        + b**2 * 2*(x[0]*sin - x[1]*cos) * (+x[0]*cos + x[1]*sin),
        ]).T


# pylint: disable=too-many-arguments
def get_ecc(f, centers, points):
    '''Get eccentricities of ellipses.
    '''
    A = []
    B = []
    for center, pts in zip(centers, points):
        X = (f(pts).T - f(center)).T

        (a, b, _), _ = leastsq(
            lambda a_b_theta: f_ellipse(a_b_theta, X),
            [1.0, 1.0, 0.0],
            Dfun=lambda a_b_theta: jac_ellipse(a_b_theta, X),
            )

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

    return numpy.array([A, B]).reshape(-1)
    # return numpy.log(ab / target_radius)


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
        center = [datak['x'], datak['y']]
        centers.append(center)
        offset = (
            numpy.array([numpy.ones(delta_y_delta_x.shape[0]), delta_y_delta_x])
            / numpy.sqrt(1 + delta_y_delta_x**2) * delta_s
            )
        points.append(numpy.column_stack([
            (center + offset.T).T,
            (center - offset.T).T,
            ]))


    centers = numpy.array(centers)

    # shift white point to center (0.0, 0.0)
    pade2d = Pade2d([2, 2, 2, 2], [1/3, 1/3])


    def f2(alpha):
        pade2d.alpha = alpha
        ecc = get_ecc(pade2d.eval, centers, points)
        # compute standard deviation of ecc
        average = numpy.sum(ecc) / len(ecc)
        out = (ecc - average) / average
        print(numpy.sum(out**2))
        return out

    # Create the identity function as initial guess
    print('num parameters: {}'.format(pade2d.total_num_coefficients))
    print('\ninitial parameters:')
    pade2d.print()
    # out = leastsq(f, coeff0, full_output=True)
    # print(out)
    # exit(1)
    # coeff1, _ = leastsq(f2, coeff0, maxfev=10000)

    ecc0 = get_ecc(pade2d.eval, centers, points)

    # Levenberg-Marquardt (lm) is better suited for small, dense, unconstrained
    # problems, but it needs more conditions than parameters.
    out = least_squares(f2, pade2d.alpha, method='trf')
    coeff1 = out.x
    print('\noptimal parameters:')
    pade2d.print()

    # plot statistics
    plt.plot(ecc0, label='ecc before')
    ecc1 = get_ecc(pade2d.eval, centers, points)
    plt.plot(ecc1, label='ecc opt')
    plt.legend()
    plt.grid()

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
        xy_to_2d=pade2d.eval,
        plot_standard_deviations=True
        )

    plt.show()
    return


if __name__ == '__main__':
    _main()
