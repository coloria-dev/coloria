import pathlib

import dufte
import matplotlib.pyplot as plt
import numpy
import yaml

from ..helpers import _compute_straight_line_residuals


def load():
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "fairchild_chen.yaml") as f:
        data = yaml.safe_load(f)

    for key in data:
        data[key]["lightness"] = numpy.asarray(data[key]["lightness"])
        data[key]["xyz"] = numpy.asarray(data[key]["xyz"])
    return data


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
    d = load()
    data = d["SL2"]

    coords = cs.from_xyz100(data["xyz"].T)
    l = coords[cs.k0]

    plt.style.use(dufte.style)

    plt.plot(l, data["lightness"], "o")
    plt.xlabel(cs.labels[cs.k0])
    plt.ylabel("experimental lightness")
    plt.title(f"Fairchild-Chen SL2 data for {cs.name}")
    plt.show()


def residuals(cs):
    wp, d = load()
    return _compute_straight_line_residuals(cs, wp, d)


def stress(cs):
    return 100 * residuals(cs)
