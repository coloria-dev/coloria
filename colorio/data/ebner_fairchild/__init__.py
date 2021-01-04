import pathlib

import matplotlib.pyplot as plt
import numpy
import yaml

from ..helpers import _compute_straight_line_residuals, _plot_color_constancy_data


def _load_data():
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "ebner_fairchild.yaml") as f:
        data = yaml.safe_load(f)

    wp = numpy.array(data["white point"])
    d = [
        numpy.column_stack([dat["reference xyz"], numpy.array(dat["same"]).T])
        for dat in data["data"]
    ]
    return wp, d


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
    wp, d = _load_data()
    _plot_color_constancy_data(d, wp, cs)


def residuals(cs):
    wp, d = _load_data()
    return _compute_straight_line_residuals(cs, wp, d)
