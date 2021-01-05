"""
David L. MacAdam,
Visual Sensitivities to Color Differences in Daylight,
Journal of the Optional Society of America,
Volume 32, May, 1942, Number 5,
https://doi.org/10.1364/JOSA.32.000247
"""
import pathlib

import matplotlib.pyplot as plt
import numpy
import yaml

from ..helpers import _plot_ellipses


def _load_data():
    # Extract ellipse centers and offsets from MacAdams data
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "table3.yaml") as f:
        data = yaml.safe_load(f)
    #
    # collect the ellipse centers and offsets
    xy_centers = []
    xy_offsets = []
    for datak in data:
        # collect ellipse points
        _, _, _, _, delta_y_delta_x, delta_s = numpy.array(datak["data"]).T
        offset = (
            numpy.array([numpy.ones(delta_y_delta_x.shape[0]), delta_y_delta_x])
            / numpy.sqrt(1 + delta_y_delta_x ** 2)
            * delta_s
        )
        if offset.shape[1] < 2:
            continue
        xy_centers.append([datak["x"], datak["y"]])
        xy_offsets.append(numpy.column_stack([+offset, -offset]))
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
    xy_centers, xy_offsets = _load_data()
    _plot_ellipses(xy_centers, xy_offsets, *args, **kwargs)


def residuals(cs):
    wp, d = _load_data()
    return _compute_straight_line_residuals(cs, wp, d)
