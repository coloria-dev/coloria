import pathlib

import matplotlib.pyplot as plt
import numpy as np
import yaml

from ...cs import CIELAB


def _load_data():
    this_dir = pathlib.Path(__file__).resolve().parent

    with open(this_dir / "berns.yaml") as f:
        data = yaml.safe_load(f)

    d = {}
    for key, value in data.items():
        value = np.array([row[2:] for row in value])
        # the vectors are already normalized, but only given in low precision. Normalize
        # them to full precision.
        vectors = value[:, 8:11]
        vectors = (vectors.T / np.linalg.norm(vectors, axis=1)).T
        d[key.lower()] = {
            "centers": value[:, 5:8],
            "vectors": vectors,
            "t50": value[:, 0],
        }
    return d


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
    d = _load_data()

    if key not in d:
        string = ", ".join(d.keys())
        raise KeyError(f"Choose one of {string}.")

    centers = d[key]["centers"]
    vectors = d[key]["vectors"]
    t50 = d[key]["t50"]

    endpoints = centers + (t50 * vectors.T).T

    cielab = CIELAB()
    cs_centers = cs.from_xyz100(cielab.to_xyz100(centers.T))
    cs_endpoints = cs.from_xyz100(cielab.to_xyz100(endpoints.T))

    # reorder the coords such that the lightness in the last (the z-)component
    cs_centers = np.roll(cs_centers, 2 - cs.k0, axis=0)
    cs_endpoints = np.roll(cs_endpoints, 2 - cs.k0, axis=0)
    labels = np.roll(cs.labels, 2 - cs.k0, axis=0)

    ax = plt.axes(projection="3d")
    # # Plot the tile points.
    # pts_2d = np.delete(pts, cs.k0, axis=0)
    for c, e in zip(cs_centers.T, cs_endpoints.T):
        rgb = cs.to_rgb1(c)
        if np.all((0 <= rgb) & (rgb <= 1)):
            color = rgb
        else:
            color = "k"
        plt.plot([c[0]], [c[1]], [c[2]], "o", color=color)
        plt.plot([c[0], e[0]], [c[1], e[1]], [c[2], e[2]], "-", color=color)

    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])


def residual(cs):
    d = _load_data()

    delta = []
    cielab = CIELAB()
    for key in d:
        centers = d[key]["centers"]
        vectors = d[key]["vectors"]
        t50 = d[key]["t50"]

        endpoints = centers + (t50 * vectors.T).T

        cs_centers = cs.from_xyz100(cielab.to_xyz100(centers.T))
        cs_endpoints = cs.from_xyz100(cielab.to_xyz100(endpoints.T))

        diff = cs_centers - cs_endpoints
        delta.append(np.linalg.norm(diff, axis=0))

    delta = np.concatenate(delta)

    diff = np.average(delta) - delta
    return np.sqrt(np.dot(diff, diff) / np.dot(delta, delta))


def stress(cs):
    return 100 * residual(cs)
