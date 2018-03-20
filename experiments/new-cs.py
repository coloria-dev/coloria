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



class MacAdam(object):
    def __init__(self):
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

        self.centers = numpy.array(centers)
        self.points = points
        return

    def get_ellipse_axes(self, f):
        '''Get eccentricities of ellipses.
        '''
        X = [
            (f.eval(pts).T - f.eval(center)).T
            for center, pts in zip(self.centers, self.points)
            ]

        ab = numpy.array([
            # Solve least squares problem for [1/a, 1/b, theta] and pick [a, b]
            1 / leastsq(
                lambda a_b_theta: f_ellipse(a_b_theta, x),
                [1.0, 1.0, 0.0],
                Dfun=lambda a_b_theta: jac_ellipse(a_b_theta, x),
                )[0][:2]
            for x in X
            ])

        return ab.flatten()
        # return numpy.log(ab / target_radius)

    def cost(self, f):
        ecc = self.get_ellipse_axes(f)
        # target = numpy.sum(ecc) / len(ecc)
        target = 0.002
        out = (ecc - target) / target
        print(numpy.sum(out**2))
        return out



class MacAdam2(object):
    def __init__(self):
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

        self.centers = numpy.array(centers)

        self.J = numpy.array(self.get_local_linearizations1(centers, points))
        # self.J = numpy.array(self.get_local_linearizations2(centers, points))

        # # plot
        # for center, pts, j in zip(centers, points, self.J):
        #     # plot points
        #     p = (pts.T - center).T
        #     plt.plot(*p, '.')
        #     # plot circle
        #     t = numpy.linspace(0.0, 2.0*numpy.pi, 1000)
        #     xy = numpy.array([numpy.cos(t), numpy.sin(t)])
        #     plt.plot(*numpy.dot(j, xy), '-', label='ellipse')
        #     plt.legend()
        #     # # plot transformation
        #     # xy_new = numpy.dot(j, p)
        #     # plt.plot(*xy_new, 'x')
        #     plt.axis('equal')
        #     plt.show()
        return

    def get_ellipse_axes(self, f):
        jacs = f.jac(self.centers.T)
        # jacs is of shape (2, 2, k); for svd, it needs shape (k, 2, 2)
        jacs = numpy.moveaxis(jacs, -1, 0)

        M = numpy.array([
            numpy.dot(jac, j) for jac, j in zip(jacs, self.J)
            ])

        # a = (M[:, 0, 0] + M[:, 1, 1]) / 2
        # b = (M[:, 0, 0] - M[:, 1, 1]) / 2
        # c = (M[:, 1, 0] + M[:, 0, 1]) / 2
        # d = (M[:, 1, 0] - M[:, 0, 1]) / 2
        # q = numpy.sqrt(a**2 + d**2)
        # r = numpy.sqrt(b**2 + c**2)
        # sigma = numpy.array([q+r, q-r]).T

        _, sigma, _ = numpy.linalg.svd(M)

        # The singular values of (invJ, jacs) are the inverses of the axis
        # lengths of the ellipses after tranformation by f.
        return sigma.flatten()

    def cost(self, f):
        ax = self.get_ellipse_axes(f)
        # target = numpy.sum(ax) / len(ax)
        target = 0.002
        out = (ax - target) / target
        print(numpy.sum(out**2))
        return out

    def get_local_linearizations1(self, centers, points):
        # Get ellipse parameters
        X = [
            (pts.T - center).T
            for center, pts in zip(centers, points)
            ]
        a_b_theta = numpy.array([
            # Solve least squares problem for [1/a, 1/b, theta]
            # and pick [a, b, theta]
            leastsq(
                lambda a_b_theta: f_ellipse(a_b_theta, x),
                [1.0, 1.0, 0.0],
                Dfun=lambda a_b_theta: jac_ellipse(a_b_theta, x),
                )[0]
            for x in X
            ])
        a_b_theta = numpy.array([
            1 / a_b_theta[:, 0],
            1 / a_b_theta[:, 1],
            a_b_theta[:, 2]
            ]).T
        # Construct 2x2 matrices that approximately convert unit circles into
        # the ellipse defined by the points.
        # rad = numpy.sum(a_b_theta[:, :2]) / (2*len(a_b_theta))
        J = []
        for abt in a_b_theta:
            a, b, theta = abt
            J.append(numpy.array([
                [a * numpy.cos(theta), -b * numpy.sin(theta)],
                [a * numpy.sin(theta), b * numpy.cos(theta)],
                ]))

        return J

    def get_local_linearizations2(self, centers, points):
        X = [
            (pts.T - center).T
            for center, pts in zip(centers, points)
            ]

        def f_linear_function(j, x):
            Jx = numpy.dot(j.reshape(2, 2), x)
            out = numpy.einsum('ij,ij->j', Jx, Jx) - 1.0
            return out

        def jac_linear_function(j, x):
            J = j.reshape(2, 2)
            return numpy.array([
                2*J[0, 0]*x[0]**2 + 2*J[0, 1]*x[0]*x[1],
                2*J[0, 1]*x[1]**2 + 2*J[0, 0]*x[0]*x[1],
                2*J[1, 0]*x[0]**2 + 2*J[1, 1]*x[0]*x[1],
                2*J[1, 1]*x[1]**2 + 2*J[1, 0]*x[0]*x[1],
                ]).T

        J = []
        for x in X:
            j, _ = leastsq(
                lambda J: f_linear_function(J, x),
                [1.0, 0.0, 0.0, 1.0],
                Dfun=lambda J: jac_linear_function(J, x),
                # full_output=True
                )
            J.append(numpy.linalg.inv(j.reshape(2, 2)))

        return J


def _main():
    pade2d = Pade2d([1, 1, 1, 1])

    # xy = numpy.random.rand(2, 5)
    # print(pade2d.eval(xy).shape)
    # exit(1)

    # macadam = MacAdam()
    macadam = MacAdam2()


    def f(alpha):
        pade2d.set_alpha(alpha)
        return macadam.cost(pade2d)

    # Create the identity function as initial guess
    print('num parameters: {}'.format(pade2d.total_num_coefficients))
    print('\ninitial parameters:')
    pade2d.print()

    alpha0 = pade2d.alpha.copy()

    # Levenberg-Marquardt (lm) is better suited for small, dense, unconstrained
    # problems, but it needs more conditions than parameters. This is not the
    # case for larger polynomial degrees.
    out = least_squares(f, pade2d.alpha, method='lm')
    pade2d.set_alpha(out.x)
    print('\noptimal parameters:')
    pade2d.print()

    # plot statistics
    pade2d.set_alpha(alpha0)
    axes0 = macadam.get_ellipse_axes(pade2d)
    plt.plot(axes0, label='axes lengths before')
    pade2d.set_alpha(out.x)
    axes1 = macadam.get_ellipse_axes(pade2d)
    plt.plot(axes1, label='axes lengths opt')
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
