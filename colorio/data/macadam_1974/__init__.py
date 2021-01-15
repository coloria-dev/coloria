"""
David L. MacAdam,
Uniform color scales,
Journal of the Optical Society of America, Vol. 64, Issue 12, pp. 1691-1702,
1974,
<https://doi.org/10.1364/JOSA.64.001691>.
"""
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import yaml

from ...cs import XYY


def _load_data():
    # Extract ellipse centers and offsets from MacAdams data
    this_dir = pathlib.Path(__file__).resolve().parent

    with open(this_dir / "table2.yaml") as f:
        data = yaml.safe_load(f)

    t = dict(zip(data.keys(), range(len(data))))
    xyy1_tiles = np.array([[val[0], val[1], val[2]] for val in data.values()])
    xyz100_tiles = XYY(100).to_xyz100(xyy1_tiles.T)

    with open(this_dir / "table1.yaml") as f:
        data = yaml.safe_load(f)

    d = np.array([item[3] for item in data])
    pairs = np.array([[t[item[1]], t[item[2]]] for item in data])
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
    d = d[np.all(pairs <= 43, axis=1)]
    pairs = pairs[np.all(pairs <= 43, axis=1)]
    xyz100_tiles = xyz100_tiles[:, :43]
    pts = cs.from_xyz100(xyz100_tiles)

    # Plot the tile points.
    pts_2d = np.delete(pts, cs.k0, axis=0)
    plt.plot(pts_2d[0], pts_2d[1], "ok", fillstyle="none")

    # for k, pt in enumerate(pts.T):
    #     # plt.text(pt[0], pt[1], k + 1)
    #     plt.plot(pt[0], pt[1], "ok", fillstyle="none")

    # scale the distances
    diff = pts[:, pairs]
    delta = np.linalg.norm(diff[..., 0] - diff[..., 1], axis=0)
    alpha = np.dot(d, delta) / np.dot(d, d)
    d *= alpha

    # plot arrow
    for dist, pair in zip(d, pairs):
        # arrow from pair[0] to pair[1]
        base = pts[:, pair[0]]
        diff = pts[:, pair[1]] - pts[:, pair[0]]
        v = diff / np.linalg.norm(diff, 2) * dist / 2
        base = np.delete(base, cs.k0)
        v = np.delete(v, cs.k0)
        plt.arrow(base[0], base[1], v[0], v[1], length_includes_head=True, color="k")
        # arrow from pair[1] to pair[0]
        base = pts[:, pair[1]]
        base = np.delete(base, cs.k0)
        v = -v
        plt.arrow(base[0], base[1], v[0], v[1], length_includes_head=True, color="k")

    plt.gca().set_aspect("equal")
    plt.show()


def residual(cs):
    d, pairs, xyz100_tiles = _load_data()

    pts = cs.from_xyz100(xyz100_tiles)

    diff = pts[:, pairs]
    delta = np.linalg.norm(diff[..., 0] - diff[..., 1], axis=0)

    alpha = np.dot(d, delta) / np.dot(d, d)
    val = np.dot(alpha * d - delta, alpha * d - delta) / np.dot(delta, delta)
    return np.sqrt(val)


def stress(cs):
    return 100 * residual(cs)
