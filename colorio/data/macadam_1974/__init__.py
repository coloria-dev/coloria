"""
David L. MacAdam,
Uniform color scales,
Journal of the Optical Society of America, Vol. 64, Issue 12, pp. 1691-1702,
1974,
<https://doi.org/10.1364/JOSA.64.001691>.
"""
import pathlib

import matplotlib.pyplot as plt
import numpy
import yaml

from ...cs import XYY1


def _load_data():
    # Extract ellipse centers and offsets from MacAdams data
    this_dir = pathlib.Path(__file__).resolve().parent

    with open(this_dir / "table2.yaml") as f:
        data = yaml.safe_load(f)

    t = dict(zip(data.keys(), range(len(data))))
    xyy1_tiles = numpy.array([[val[0], val[1], val[2]] for val in data.values()])
    xyy1_tiles[:, 2] /= 100
    xyz100_tiles = XYY1().to_xyz100(xyy1_tiles.T)

    with open(this_dir / "table1.yaml") as f:
        data = yaml.safe_load(f)

    d = numpy.array([item[3] for item in data])
    pairs = numpy.array([[t[item[1]], t[item[2]]] for item in data])
    return d, pairs, xyz100_tiles


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


def plot(cs):
    d, pairs, xyz100_tiles = _load_data()

    # Plot the tiles points. Scratch the lightness component.
    xyz100_tiles = xyz100_tiles[:, :43]
    print(xyz100_tiles.shape)
    pts = cs.from_xyz100(xyz100_tiles)
    pts = numpy.delete(pts, cs.k0, axis=0)

    for k, pt in enumerate(pts.T):
        plt.text(pt[0], pt[1], k + 1)
        plt.plot(pt[0], pt[1], "ok", fillstyle="none")
    # plt.plot(pts[0], pts[1], "ok", fillstyle="none")
    plt.gca().set_aspect("equal")
    plt.show()

    # plot arrows
    for dist, pair in zip(d, pairs):
        print(pair)
        exit(1)

    # plt.plt(pts[0
    # xy_centers, xy_offsets = _load_data()
    # _plot_ellipses(xy_centers, xy_offsets, *args, ellipse_scaling, **kwargs)


def residual(cs):
    d, pairs, xyz100_tiles = _load_data()

    pts = cs.from_xyz100(xyz100_tiles)

    diff = pts[:, pairs]
    delta = numpy.linalg.norm(diff[..., 0] - diff[..., 1], axis=0)

    alpha = numpy.dot(d, delta) / numpy.dot(d, d)
    val = numpy.dot(alpha * d - delta, alpha * d - delta) / numpy.dot(delta, delta)
    return numpy.sqrt(val)
