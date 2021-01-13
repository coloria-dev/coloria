import pathlib

import matplotlib.pyplot as plt
import numpy
import yaml

from ...illuminants import whitepoints_cie1931
from ..helpers import _compute_straight_line_residuals, _plot_hue_linearity_data


def load():
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "table3.yaml") as f:
        data = yaml.safe_load(f)

    wp = whitepoints_cie1931["C"]
    d = [numpy.array(list(color.values())).T for color in data.values()]
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
    wp, d = load()
    _plot_hue_linearity_data(d, wp, cs)


def residuals(cs):
    wp, d = load()
    return _compute_straight_line_residuals(cs, wp, d)


def stress(cs):
    return 100 * residuals(cs)
