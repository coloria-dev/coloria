"""
M. R. Luo, B. Rigg,
Chromaticity Discrimination Ellipses for Surface Colours,
Color Research and Application, Volume 11, Issue 1, Spring 1986, Pages 25-42,
<https://doi.org/10.1002/col.5080110107>.
"""
import pathlib

import matplotlib.pyplot as plt
import numpy
import yaml

from ..helpers import _compute_ellipse_residual, _plot_ellipses


def load(num_offset_points: int):
    # Extract ellipse centers and offsets from MacAdams data
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "luo-rigg.yaml") as f:
        data = yaml.safe_load(f)
    #
    xyy100_centers = []
    xyy100_points = []
    # collect the ellipse centers and offsets
    alpha = numpy.linspace(0.0, 2 * numpy.pi, num_offset_points, endpoint=False)
    circle_pts = numpy.array([numpy.cos(alpha), numpy.sin(alpha)])
    for data_set in data.values():
        # The set factor is the mean of the R values
        # set_factor = sum([dat[-1] for dat in data_set.values()]) / len(data_set)
        for x, y, Y, a, a_div_b, theta_deg, _ in data_set.values():
            theta = theta_deg * 2 * numpy.pi / 360
            a /= 1.0e4
            a *= (Y / 30) ** 0.2
            b = a / a_div_b
            # plot the ellipse
            c = numpy.array([x, y])
            J = numpy.array(
                [
                    [+a * numpy.cos(theta), -b * numpy.sin(theta)],
                    [+a * numpy.sin(theta), +b * numpy.cos(theta)],
                ]
            )
            offsets = numpy.dot(J, circle_pts)
            pts = (c + offsets.T).T

            xyy100_centers.append(numpy.array([*c, Y]))
            xyy100_points.append(numpy.array([*pts, numpy.full(pts.shape[1], Y)]))

    return xyy100_centers, xyy100_points


def show(*args, **kwargs):
    plt.figure()
    plot(*args, **kwargs)
    plt.show()
    plt.close()


def savefig(filename, *args, **kwargs):
    plt.figure()
    plot(*args, **kwargs)
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


def plot(cs, ellipse_scaling: float = 2.0, num_offset_points: int = 16):
    xyy100_centers, xyy100_points = load(num_offset_points)
    _plot_ellipses(xyy100_centers, xyy100_points, cs, ellipse_scaling)
    plt.title(f"Luo-Rigg ellipses for {cs.name}")


def residuals(cs, num_offset_points):
    xyy100_centers, xyy100_points = load(num_offset_points)
    return _compute_ellipse_residual(cs, xyy100_centers, xyy100_points)


def stress(*args):
    return 100 * residuals(*args)
