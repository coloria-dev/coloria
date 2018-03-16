# -*- coding: utf-8 -*-
#
import os

import colorio
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy
from scipy.optimize import leastsq, least_squares
import yaml


class New(object):
    def __init__(self):
        a = colorio.illuminants.whitepoints_cie1931['D65']
        b = numpy.array([100.0, 0.0, 0.0])
        # Get a matrix that maps a to b.
        # <https://math.stackexchange.com/a/2672702/36678>
        aa = a / numpy.sqrt(numpy.dot(a, a))
        bb = b / numpy.sqrt(numpy.dot(b, b))
        ab = aa + bb
        self.M = 2 * numpy.outer(ab, ab) / numpy.dot(ab, ab) - numpy.eye(3)
        self.M *= numpy.sqrt(numpy.dot(b, b)) / numpy.sqrt(numpy.dot(a, a))
        return

    def from_xyz100(self, xyz):
        return numpy.dot(self.M, xyz)


def _main():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, '../colorio/data/macadam1942/table3.yaml')) as f:
        data = yaml.safe_load(f)

    ax = plt.gca()

    scaling = 10

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


    # # plot the points
    # for center, pts in zip(centers, points):
    #     plt.plot(*center, 'x', color='k')
    #     plt.plot(*pts, '.', color='k')
    # plt.gca().set_aspect('equal')
    # plt.show()

    centers = numpy.array(centers)

    target_radius = 2.0e-3


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


    def transform(xy, alpha):
        alpha1 = numpy.concatenate([[0.0], alpha[:27]])
        a1 = create_triangle(alpha1, 6)

        alpha2 = numpy.concatenate([[1.0], alpha[27:54]])
        a2 = create_triangle(alpha2, 6)

        beta1 = numpy.concatenate([[0.0], alpha[54:81]])
        b1 = create_triangle(beta1, 6)

        beta2 = numpy.concatenate([[1.0], alpha[81:]])
        b2 = create_triangle(beta2, 6)

        return numpy.array([
            evaluate_2d_polynomial(xy, a1) / evaluate_2d_polynomial(xy, a2),
            evaluate_2d_polynomial(xy, b1) / evaluate_2d_polynomial(xy, b2),
            ])


    def f1(alpha):
        '''Function that returns the difference between target_radius and the
        distances of the ellipse points to the centers. Minimizing this
        function should make the MacAdams ellipses circles. However, what may
        also happen is that the ellipse points are lined up at target_radius
        away from the centers, leading to very sharp ellipses instead.
        '''
        dist = []
        for center, pts in zip(centers, points):
            tcenter = transform(center, alpha)
            tpts = transform(pts, alpha)

            for pt in tpts.T:
                diff = tcenter - pt
                dist.append(numpy.sqrt(numpy.dot(diff, diff)))

        dist = numpy.array(dist)
        # plt.plot([0, len(dist)], [target_radius, target_radius])
        # plt.plot(numpy.arange(len(dist)), dist)
        # plt.show()
        return abs(dist - target_radius)

    def get_radii(alpha):
        A = []
        B = []
        for center, pts in zip(centers, points):
            tcenter = transform(center, alpha)
            X = (transform(pts, alpha).T - tcenter).T

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

            (a, b, theta), _ = \
                leastsq(f_ellipse, [1.0, 1.0, 0.0], Dfun=jac_ellipse)

            A.append(1/a)
            B.append(1/b)

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

    def f2(alpha):
        return get_radii(alpha) - target_radius


    coeff0 = numpy.zeros(108)
    coeff0[0] = 1.0
    coeff0[55] = 1.0
    print(coeff0)
    # out = leastsq(f, coeff0, full_output=True)
    # print(out)
    # exit(1)
    # coeff1, _ = leastsq(f2, coeff0, maxfev=10000)

    # Levenberg-Marquardt (lm) is better suited for small, dense, unconstrained
    # problems, but it needs more conditions than parameters.
    out = least_squares(f2, coeff0, method='trf')
    coeff1 = out.x
    print(coeff1)

    radii0 = get_radii(coeff0)
    radii1 = get_radii(coeff1)
    plt.plot(numpy.arange(len(radii0)), radii0, label='radii0')
    plt.plot(numpy.arange(len(radii1)), radii1, label='radii1')
    plt.legend()
    plt.show()

    plt.figure()
    colorio.plot_macadam(
        scaling=10,
        plot_standard_deviations=True
        )
    plt.figure()
    colorio.plot_macadam(
        scaling=10,
        xy_to_2d=lambda xy: transform(xy, coeff1),
        plot_standard_deviations=True
        )
    plt.show()
    return


if __name__ == '__main__':
    _main()
