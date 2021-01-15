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

from ..._exceptions import ColorioError
from ...cs import XYY


def _load_data(remove_yellow_8=True):
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

    if remove_yellow_8:
        xyy_samples["yellow"][8] = numpy.nan

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

    if remove_yellow_8:
        distances["yellow"][numpy.any(pairs == 8, axis=1)] = numpy.nan

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
    # only plot one tile set for now
    xyy_samples, _, _ = _load_data()

    if key not in xyy_samples:
        string = ", ".join(xyy_samples.keys())
        raise ColorioError(f"`key` must be one of {string}.")

    xyy = xyy_samples[key]
    # kick out the nans
    xyy = xyy[~numpy.any(numpy.isnan(xyy), axis=1)]
    xyz100 = XYY(100).to_xyz100(xyy.T)
    coords = cs.from_xyz100(xyz100)

    # reorder the coords such that the lightness in the last (the z-)component
    coords = numpy.roll(coords, 2 - cs.k0, axis=0)
    labels = numpy.roll(cs.labels, 2 - cs.k0, axis=0)

    ax = plt.axes(projection="3d")
    ax.scatter(*coords, marker="o", color=key)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])
    ax.set_title(f"Witt dataset, {key} tiles, in {cs.name}")


def residual(cs):
    xyy_samples, pairs, target_distances = _load_data()

    delta = []
    d = []
    for key in xyy_samples:
        isnan = numpy.isnan(target_distances[key])
        d.append(target_distances[key][~isnan])
        # compute the actual distances in the color space `cs`
        xyz_samples = XYY(100).to_xyz100(xyy_samples[key].T)
        cs_samples = cs.from_xyz100(xyz_samples).T
        cs_diff = cs_samples[pairs[~isnan, 0]] - cs_samples[pairs[~isnan, 1]]
        cs_dist = numpy.sqrt(numpy.einsum("ij,ij->i", cs_diff, cs_diff))
        delta.append(cs_dist)

    d = numpy.concatenate(d)
    delta = numpy.concatenate(delta)

    assert d.shape == delta.shape
    # The original article lists 418 pairs, but Yellow-8 is most likely printed with an
    # error in the article. After its removal, 414 sample pairs remain.
    assert len(d) == 414

    alpha = numpy.dot(d, delta) / numpy.dot(d, d)
    val = numpy.dot(alpha * d - delta, alpha * d - delta) / numpy.dot(delta, delta)
    return numpy.sqrt(val)


def stress(cs):
    return 100 * residual(cs)
