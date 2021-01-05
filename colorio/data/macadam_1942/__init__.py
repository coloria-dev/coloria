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

from ..._helpers import _plot_ellipses


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


def plot(*args, **kwargs):
    xy_centers, xy_offsets = _load_data()
    _plot_ellipses(xy_centers, xy_offsets, *args, **kwargs)


def residuals(cs, Y: float):
    distances = []
    xy_centers, xy_offsets = _load_data()

    # scale invariance by normalizing on center average
    xy_centers = numpy.array(xy_centers).T
    xyy_centers = numpy.array([*xy_centers, numpy.full(xy_centers.shape[1], Y)])
    xyz_centers = _xyy_to_xyz100(xyy_centers)
    cs_centers = cs.from_xyz100(xyz_centers)
    cs_centers = numpy.delete(cs_centers, cs.k0, axis=0)
    avg = numpy.average(cs_centers, axis=1)
    alpha = numpy.linalg.norm(avg)

    for xy_center, xy_offsets in zip(xy_centers.T, xy_offsets):
        xy_ellips = (xy_center + xy_offsets.T).T
        # append Y
        xyy_center = numpy.array([*xy_center, Y])
        xyy_ellips = numpy.array([*xy_ellips, numpy.full(xy_ellips.shape[1], Y)])
        xyz_center = _xyy_to_xyz100(xyy_center)
        assert numpy.all(xyz_center >= 0.0)
        xyz_ellips = _xyy_to_xyz100(xyy_ellips)
        assert numpy.all(xyz_ellips >= 0.0)
        # plt.plot(xyz_center[0], xyz_center[2], "x")
        # plt.plot(xyz_ellips[0], xyz_ellips[2], "o")
        # plt.show()
        cs_center = cs.from_xyz100(xyz_center)
        cs_ellips = cs.from_xyz100(xyz_ellips)
        # remove lightness data
        cs_center = numpy.delete(cs_center, cs.k0, axis=0)
        cs_ellips = numpy.delete(cs_ellips, cs.k0, axis=0)
        # plt.plot(cs_center[0], cs_center[1], "x")
        # plt.plot(cs_ellips[0], cs_ellips[1], "o")
        # plt.show()
        #
        # scale and compute distances to center
        cs_center /= alpha
        cs_ellips /= alpha
        diff = (cs_center - cs_ellips.T).T
        distances.append(numpy.sqrt(diff[0] ** 2 + diff[1] ** 2))

    distances = numpy.concatenate(distances)
    avg = numpy.average(distances)
    return numpy.sqrt(numpy.sum((distances - avg) ** 2))


def _xyy_to_xyz100(xyy):
    x, y, Y = xyy
    return numpy.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100
