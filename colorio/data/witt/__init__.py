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
    this_dir = pathlib.Path(__file__).resolve().parent

    with open(this_dir / "table_a2.yaml") as f:
        data = yaml.safe_load(f)

    print(data)
    exit(1)

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

    # only consider the first 43 tiles which are all of approximately the same lightness
    # TODO plot 3D tetrahedral pairs, too
    d = d[numpy.all(pairs <= 43, axis=1)]
    pairs = pairs[numpy.all(pairs <= 43, axis=1)]
    xyz100_tiles = xyz100_tiles[:, :43]
    pts = cs.from_xyz100(xyz100_tiles)

    # Plot the tile points.
    pts_2d = numpy.delete(pts, cs.k0, axis=0)
    plt.plot(pts_2d[0], pts_2d[1], "ok", fillstyle="none")

    # for k, pt in enumerate(pts.T):
    #     # plt.text(pt[0], pt[1], k + 1)
    #     plt.plot(pt[0], pt[1], "ok", fillstyle="none")

    # scale the distances
    diff = pts[:, pairs]
    delta = numpy.linalg.norm(diff[..., 0] - diff[..., 1], axis=0)
    alpha = numpy.dot(d, delta) / numpy.dot(d, d)
    d *= alpha

    # plot arrow
    for dist, pair in zip(d, pairs):
        # arrow from pair[0] to pair[1]
        base = pts[:, pair[0]]
        diff = pts[:, pair[1]] - pts[:, pair[0]]
        v = diff / numpy.linalg.norm(diff, 2) * dist / 2
        base = numpy.delete(base, cs.k0)
        v = numpy.delete(v, cs.k0)
        plt.arrow(base[0], base[1], v[0], v[1], length_includes_head=True, color="k")
        # arrow from pair[1] to pair[0]
        base = pts[:, pair[1]]
        base = numpy.delete(base, cs.k0)
        v = -v
        plt.arrow(base[0], base[1], v[0], v[1], length_includes_head=True, color="k")

    plt.gca().set_aspect("equal")
    plt.show()


def residual(cs):
    d, pairs, xyz100_tiles = _load_data()

    pts = cs.from_xyz100(xyz100_tiles)

    diff = pts[:, pairs]
    delta = numpy.linalg.norm(diff[..., 0] - diff[..., 1], axis=0)

    alpha = numpy.dot(d, delta) / numpy.dot(d, d)
    val = numpy.dot(alpha * d - delta, alpha * d - delta) / numpy.dot(delta, delta)
    return numpy.sqrt(val)


def stress(cs):
    return 100 * residual(cs)
