"""
Klaus Witt,
Geometric relations between scales of small colour differences,
Color Research and Application, Volume 24, Issue 2, April 1999, Pages 78-92,
<https://doi.org/10.1002/(SICI)1520-6378(199904)24:2<78::AID-COL3>3.0.CO;2-M>.
"""
import pathlib

import matplotlib.pyplot as plt
import numpy
import yaml

from ...cs import XYY1


def _load_data():
    this_dir = pathlib.Path(__file__).resolve().parent

    with open(this_dir / "table_a1.yaml") as f:
        xyy_samples = yaml.safe_load(f)
    xyy_samples = numpy.array([s for s in xyy_samples.values()])
    xyy_samples = {
            "yellow": xyy_samples[:, 0],
            "grey": xyy_samples[:, 1],
            "green": xyy_samples[:, 2],
            "red": xyy_samples[:, 3],
            "blue": xyy_samples[:, 4],
    }

    with open(this_dir / "table_a2.yaml") as f:
        data = yaml.safe_load(f)

    # each line has 12 entries:
    # pair, yellow (mean + sigma), grey (m+s), green (m+s), red (m+s), blue (m+s)
    pairs = numpy.array([item[:2] for item in data])
    distances = {
        "yellow": numpy.array([item[2] for item in data]),
        "grey": numpy.array([item[4] for item in data]),
        "green": numpy.array([item[6] for item in data]),
        "red": numpy.array([item[8] for item in data]),
        "blue": numpy.array([item[10] for item in data]),
    }
    return xyy_samples, pairs, distances


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


def plot(cs, key):
    xyy_samples, _, _ = _load_data()

    # only plot yellow for now
    # d = distances[key]
    xyy = xyy_samples[key]

    xyz100 = XYY1().to_xyz100(xyy

    print(xyy.shape)
    exit(1)
    ax = plt.axes(projection="3d")
    ax.scatter(*xyy.T, marker="o", color=key)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("Y")
    # plt.zlabel("Y")
    # plt.title(
    # ax.set_box_aspect((1, 1, 1))


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
