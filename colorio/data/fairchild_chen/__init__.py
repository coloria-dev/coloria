import pathlib

import dufte
import matplotlib.pyplot as plt
import numpy
import yaml


def load():
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "fairchild_chen.yaml") as f:
        data = yaml.safe_load(f)

    for key in data:
        data[key]["lightness"] = numpy.asarray(data[key]["lightness"])
        data[key]["xyz"] = numpy.asarray(data[key]["xyz"])
    return data


def show(*args):
    plt.figure()
    plot(*args)
    plt.show()
    plt.close()


def savefig(filename, *args):
    plt.figure()
    plot(*args)
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


def plot(cs, key: str):
    assert key in ["SL1", "SL2"]

    data = load()[key]

    # prediced lightness
    L_ = cs.from_xyz100(data["xyz"].T)[cs.k0]

    plt.style.use(dufte.style)

    plt.plot(L_, data["lightness"], "o")
    plt.xlabel(cs.labels[cs.k0])
    plt.ylabel("experimental lightness")
    plt.title(f"Fairchild-Chen {key} data for {cs.name}")
    plt.show()


def residual(cs, key: str):
    d = load()[key]

    # experimental lightness
    L = d["lightness"]
    # predicted lightness
    L_ = cs.from_xyz100(d["xyz"].T)[cs.k0]

    alpha = numpy.dot(L, L_) / numpy.dot(L, L)
    diff = alpha * L - L_
    return numpy.sqrt(numpy.dot(diff, diff) / numpy.dot(L_, L_))


def stress(*args):
    return 100 * residual(*args)
