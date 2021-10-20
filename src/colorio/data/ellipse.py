from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike

from ..cs import XYY, ColorCoordinates, ColorSpace


class EllipseDataset:
    def __init__(self, name: str, xyy100_centers: ArrayLike, xyy100_points: ArrayLike):
        self.name = name
        self.xyy100_centers = ColorCoordinates(xyy100_centers, XYY(100))
        self.xyy100_points = ColorCoordinates(xyy100_points, XYY(100))

    def stress(self, cs: ColorSpace):
        cs_centers = self.xyy100_centers.convert(cs)
        cs_points = self.xyy100_points.convert(cs)

        diff = (cs_centers.data[:, None] - cs_points.data).reshape(3, -1)
        distances = np.sqrt(np.einsum("ij,ij->j", diff, diff))

        alpha = np.average(distances)
        return 100 * np.sqrt(np.sum((alpha - distances) ** 2) / np.sum(distances ** 2))

    def plot(self, cs: ColorSpace, ellipse_scaling: float = 1.0):
        _plot_ellipses(
            cs, self.xyy100_centers, self.xyy100_points, ellipse_scaling=ellipse_scaling
        )
        plt.title(f"{self.name} ellipses for {cs.name}")
        return plt


def _plot_ellipses(
    cs: ColorSpace,
    xyy100_centers: ColorCoordinates,
    xyy100_points: ColorCoordinates,
    ellipse_scaling: float,
):
    from matplotlib.patches import Ellipse
    from scipy.optimize import leastsq

    # make the ellipses the same color as the axes labels
    color = plt.gca().xaxis.label.get_color()

    cs_centers = xyy100_centers.convert(cs).data_hue.T
    cs_points = xyy100_points.convert(cs).data_hue.T

    for tcenter, tvals in zip(cs_centers, cs_points):
        # Given these new transformed vals, find the ellipse that best fits those
        # points
        X = (tvals - tcenter).T

        def f_ellipse(a_b_theta):
            a, b, theta = a_b_theta
            sin_t = np.sin(theta)
            cos_t = np.cos(theta)
            return (
                +(a ** 2) * (X[0] * cos_t + X[1] * sin_t) ** 2
                + b ** 2 * (X[0] * sin_t - X[1] * cos_t) ** 2
                - 1.0
            )

        def jac(a_b_theta):
            a, b, theta = a_b_theta
            x0sin = X[0] * np.sin(theta)
            x0cos = X[0] * np.cos(theta)
            x1sin = X[1] * np.sin(theta)
            x1cos = X[1] * np.cos(theta)
            return np.array(
                [
                    +2 * a * (x0cos + x1sin) ** 2,
                    +2 * b * (x0sin - x1cos) ** 2,
                    +(a ** 2) * 2 * (x0cos + x1sin) * (-x0sin + x1cos)
                    + b ** 2 * 2 * (x0sin - x1cos) * (x0cos + x1sin),
                ]
            ).T

        # We need to use some optimization here to find the new ellipses which best
        # fit the modified data.
        (a, b, theta), _ = leastsq(f_ellipse, [1.0, 1.0, 0.0], Dfun=jac)

        # plot the scaled ellipse
        e = Ellipse(
            xy=tcenter,
            width=ellipse_scaling * 2 / a,
            height=ellipse_scaling * 2 / b,
            angle=theta / np.pi * 180,
            # label=label,
        )
        plt.gca().add_patch(e)
        # e.set_alpha(0.5)

        e.set_facecolor(color)

        # plt.plot(*tcenter, "xk")
        # plt.plot(*tvals, "ok")
        # plt.show()

    plt.gca().set_aspect("equal")
    plt.xlabel(cs.hue_labels[0])
    plt.ylabel(cs.hue_labels[1], rotation=0)

    # mpl doesn't update axis limits when adding artists,
    # <https://github.com/matplotlib/matplotlib/issues/19290>.
    # Handle it manually for now.
    tcenters = xyy100_centers.convert(cs).data_hue

    xmin = np.min(tcenters[0])
    xmax = np.max(tcenters[0])
    ymin = np.min(tcenters[1])
    ymax = np.max(tcenters[1])
    width = xmax - xmin
    height = ymax - ymin
    plt.xlim(xmin - 0.2 * width, xmax + 0.2 * width)
    plt.ylim(ymin - 0.2 * height, ymax + 0.2 * height)
