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
    a -
