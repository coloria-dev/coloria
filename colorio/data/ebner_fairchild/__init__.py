import pathlib

import matplotlib.pyplot as plt
import numpy
import yaml

from ..helpers import _plot_color_constancy_data


def show(cs):
    plt.figure()
    plot(cs)
    plt.show()
    plt.close()


def savefig(cs, filename):
    plt.figure()
    plot(cs)
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


def plot(cs):
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "ebner_fairchild.yaml") as f:
        data = yaml.safe_load(f)

    wp = numpy.array(data["white point"])

    d = [
        numpy.column_stack([dat["reference xyz"], numpy.array(dat["same"]).T])
        for dat in data["data"]
    ]

    _plot_color_constancy_data(d, wp, cs)
    plt.grid(False)
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)


def residuals(cs):
    """Compute the TLS residuals for each of the arms."""
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "ebner_fairchild.yaml") as f:
        data = yaml.safe_load(f)

    wp = numpy.array(data["white point"])
    d = [
        numpy.column_stack([dat["reference xyz"], numpy.array(dat["same"]).T])
        for dat in data["data"]
    ]

    # remove the row corresponding to lightness
    idx = [0, 1, 2]
    idx.pop(cs.k0)
    wp_cs = cs.from_xyz100(wp)[idx]
    s2 = []
    for dd in d:
        vals = cs.from_xyz100(dd)[idx]
        # move values such that whitepoint is in the origin
        vals = (vals.T - wp_cs).T
        # scale by average to achieve scale invariance
        avg = numpy.sum(vals, axis=1) / vals.shape[1]
        vals /= numpy.linalg.norm(avg)
        # could also be computed explicitly
        s2.append(numpy.linalg.svd(vals, compute_uv=False)[-1])
        # plt.plot(vals[0], vals[1], "x")
        # plt.gca().set_aspect("equal")
        # plt.show()

    return s2
