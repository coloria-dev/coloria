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

    target_dist = 1.0e-3

    def transform(xy, alpha):
        x, y = xy
        return numpy.array([
            alpha[0]*x + alpha[1]*y,
            alpha[2]*x + alpha[3]*y,
            ])

    def f(alpha):
        dist = []
        for center, pts in zip(centers, points):
            tcenter = transform(center, alpha)
            tpts = transform(pts, alpha)

            for pt in tpts.T:
                diff = tcenter - pt
                dist.append(numpy.sqrt(numpy.dot(diff, diff)))

        dist = numpy.array(dist)
        # plt.plot(numpy.arange(len(dist)), dist)
        # plt.show()
        return (dist - target_dist)**2


    coeff0 = [1.0, 0.0, 0.0, 1.0]
    print(coeff0)
    coeff1, _ = leastsq(f, coeff0)
    print(coeff1)

    # # plot the points
    # for center, pts in zip(centers, points):
    #     center = transform(center, coeff0)
    #     pts = transform(pts, coeff1)
    #     plt.plot(*center, 'x', color='k')
    #     plt.plot(*pts, '.', color='k')

    colorio.show_macadam(
        scaling=10,
        xy_to_2d=lambda xy: transform(xy, coeff1),
        plot_standard_deviations=True
        )
    return


if __name__ == '__main__':
    _main()
