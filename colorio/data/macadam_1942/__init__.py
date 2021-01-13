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

from ...cs import XYY1
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
        xy_centers.append(numpy.array([datak["x"], datak["y"]]))
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


def plot(cs, ellipse_scaling=10.0):
    xy_centers, xy_offsets = _load_data()

    Y = 50.0
    xyy100_centers = []
    xyy100_points = []
    for c, off in zip(xy_centers, xy_offsets):
        xyy100_centers.append(numpy.array([*c, Y]))
        p = (c + off.T).T
        xyy100_points.append(numpy.array([*p, numpy.full(p.shape[1], Y)]))

    _plot_ellipses(xyy100_centers, xyy100_points, cs, ellipse_scaling)
    plt.title(f"MacAdam ellipses for {cs.name}")

    # cs.plot_visible_slice(
    #     lightness,
    #     outline_prec=outline_prec,
    #     fill_color=visible_gamut_fill_color,
    # )
    # if plot_srgb_gamut:
    #     cs.plot_rgb_slice(lightness)


def residual(cs, Y100: float) -> float:
    xy_centers, xy_offsets = _load_data()

    xyy1 = XYY1()

    dists = []
    for c, off in zip(xy_centers, xy_offsets):
        pts = (c + off.T).T

        # get ellipse center in transformed space
        xyY1_center = numpy.append(c, Y100 / 100)
        c = cs.from_xyz100(xyy1.to_xyz100(xyY1_center))

        # get ellipse points in transformed space
        xyY1_pts = numpy.array([*pts, numpy.full(pts.shape[1], Y100 / 100)])
        pts = cs.from_xyz100(xyy1.to_xyz100(xyY1_pts))

        # compute the distance in the transformed space
        diff = (c - pts.T).T
        dists.append(numpy.sqrt(numpy.einsum("ij,ij->j", diff, diff)))

    delta = numpy.concatenate(dists)

    alpha = numpy.average(delta)
    val = numpy.sqrt(numpy.dot(alpha - delta, alpha - delta) / numpy.dot(delta, delta))
    return val


def stress(cs, Y100: float):
    return 100 * residual(cs, Y100)
