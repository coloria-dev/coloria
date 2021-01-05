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

from ..._helpers import _plot_ellipses
from ..helpers import _compute_ellipse_residual


def _load_data(num_offset_points):
    # Extract ellipse centers and offsets from MacAdams data
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "luo-rigg.yaml") as f:
        data = yaml.safe_load(f)
    #
    xy_centers = []
    xy_offsets = []
    # collect the ellipse centers and offsets
    alpha = (
        2
        * numpy.pi
        * numpy.linspace(0.0, 2 * numpy.pi, num_offset_points, endpoint=False)
    )
    pts = numpy.array([numpy.cos(alpha), numpy.sin(alpha)])
    for data_set in data.values():
        # The set factor is the mean of the R values
        # set_factor = sum([dat[-1] for dat in data_set.values()]) / len(data_set)
        for x, y, Y, a, a_div_b, theta_deg, _ in data_set.values():
            theta = theta_deg * 2 * numpy.pi / 360
            a /= 1.0e4
            a *= (Y / 30) ** 0.2
            b = a / a_div_b
            # plot the ellipse
            xy_centers.append([x, y])
            J = numpy.array(
                [
                    [+a * numpy.cos(theta), -b * numpy.sin(theta)],
                    [+a * numpy.sin(theta), +b * numpy.cos(theta)],
                ]
            )
            xy_offsets.append(numpy.dot(J, pts))
    return xy_centers, xy_offsets


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


def plot(*args, **kwargs):
    xy_centers, xy_offsets = _load_data(num_offset_points=8)
    _plot_ellipses(xy_centers, xy_offsets, *args, **kwargs)


def residuals(cs, Y: float):
    xy_centers, xy_offsets = _load_data(num_offset_points=8)
    _compute_ellipse_residual(cs, xy_centers, xy_offsets, Y)
