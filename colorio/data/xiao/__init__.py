import pathlib

import matplotlib.pyplot as plt
import numpy as np
import yaml

from ..helpers import _compute_straight_line_residuals, _plot_hue_linearity_data


def load():
    this_dir = pathlib.Path(__file__).resolve().parent

    data = []
    with open(this_dir / "averages.yaml") as f:
        data = yaml.safe_load(f)
    data = np.array(list(data.values()))

    # Use Xiao's 'neutral gray' as white point.
    with open(this_dir / "neutral_gray.yaml") as f:
        ng_data = np.array(yaml.safe_load(f))
    ng = np.sum(ng_data, axis=0) / np.prod(ng_data.shape[:1])

    data = np.moveaxis(data, 1, 2)
    return ng, data


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


def plot(self):
    wp, d = load()
    _plot_hue_linearity_data(d, wp, self)


def residuals(cs):
    wp, d = load()
    return _compute_straight_line_residuals(cs, wp, d)


def stress(cs):
    return 100 * residuals(cs)
