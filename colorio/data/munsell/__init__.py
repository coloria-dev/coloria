import pathlib

import matplotlib.pyplot as plt
import numpy
import yaml

from ...cs import SrgbLinear


def load():
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "real.yaml") as f:
        data = yaml.safe_load(f)

    h = numpy.array(data["h"])
    V = numpy.array(data["V"])
    C = numpy.array(data["C"])
    xyy = numpy.array([data["x"], data["y"], data["Y"]])

    return h, V, C, xyy


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


def plot(cs, V):
    _, v, _, xyy = load()

    # pick the data from the given munsell level
    xyy = xyy[:, v == V]

    x, y, Y = xyy
    xyz100 = numpy.array([Y / y * x, Y, Y / y * (1 - x - y)])
    vals = cs.from_xyz100(xyz100)

    srgb = SrgbLinear()
    rgb = srgb.from_xyz100(xyz100)
    is_legal_srgb = numpy.all((0 <= rgb) & (rgb <= 1), axis=0)

    idx = [0, 1, 2]
    k1, k2 = idx[: cs.k0] + idx[cs.k0 + 1 :]

    # plot the ones that cannot be represented in SRGB in black
    plt.plot(
        vals[k1, ~is_legal_srgb],
        vals[k2, ~is_legal_srgb],
        "o",
        color="white",
        markeredgecolor="black",
    )
    # plot the srgb dots in color
    for val, rgb_ in zip(vals[:, is_legal_srgb].T, rgb[:, is_legal_srgb].T):
        plt.plot(val[k1], val[k2], "o", color=srgb.to_rgb1(rgb_))

    plt.grid()
    plt.title(f"V={V}")
    plt.xlabel(cs.labels[k1])
    plt.ylabel(cs.labels[k2])
    plt.axis("equal")
