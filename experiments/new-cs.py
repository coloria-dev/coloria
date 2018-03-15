# -*- coding: utf-8 -*-
#
import os

import colorio
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy
from scipy.optimize import leastsq
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

    def transform(xy, alpha):
        x, y = xy
        return numpy.array([
            alpha[0]*x + alpha[2]*y
            + alpha[4]*x**2 + alpha[6]*x*y + alpha[8]*y**2,
            # + alpha[10]*x**3 + alpha[12]*x**2*y + alpha[14]*x*y**2 + alpha[16]*y**3,
            alpha[1]*x + alpha[3]*y
            + alpha[5]*x**2 + alpha[7]*x*y + alpha[9]*y**2,
            # + alpha[11]*x**3 + alpha[13]*x**2*y + alpha[15]*x*y**2 + alpha[17]*y**3,
            ])

    def f1(alpha):
        '''Function that returns the difference between target_radius and the
        distances of the ellipse points to the centers. Minimizing this
        function should make the MacAdams ellipses circles. However, what may
        also happen is that the ellipse points are lined up at target_radius away
        from the centers, leading to very sharp ellipses instead.
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

    def f2(alpha):
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

        ab = numpy.concatenate([A, B])
        return abs(ab - target_radius)


    coeff0 = numpy.zeros(10)
    coeff0[0] = 1.0
    coeff0[3] = 1.0
    print(coeff0)
    # out = leastsq(f, coeff0, full_output=True)
    # print(out)
    # exit(1)
    coeff1, _ = leastsq(f2, coeff0, maxfev=10000)
    print(coeff1)

    vals0 = f2(coeff0)
    vals1 = f2(coeff1)
    plt.plot(numpy.arange(len(vals0)), vals0, label='vals0')
    plt.plot(numpy.arange(len(vals1)), vals1, label='vals1')
    plt.legend()
    plt.show()

    colorio.show_macadam(
        scaling=10,
        xy_to_2d=lambda xy: transform(xy, coeff1),
        plot_standard_deviations=True
        )
    return


if __name__ == '__main__':
    _main()
