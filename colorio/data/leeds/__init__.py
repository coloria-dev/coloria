import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np

from ..._exceptions import ColorioError
from ...cs import XYY
from ..helpers import Dataset


class Leeds(Dataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "leeds.json") as f:
            data = json.load(f)

        self.dist = np.asarray(data["dv"])
        self.xyz_pairs = np.asarray(data["pairs"])

    def plot(self, cs):
        coords = cs.from_xyz100(self.xyz_pairs.T).T

        # reorder the coords such that the lightness in the last (the z-)component
        coords = np.roll(coords, 2 - cs.k0, axis=0)
        labels = np.roll(cs.labels, 2 - cs.k0, axis=0)

        ax = plt.axes(projection="3d")
        for pair in coords:
            ax.plot(pair[:, 0], pair[:, 1], pair[:, 2], "-k")

        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_zlabel(labels[2])
        ax.set_title(f"Leeds dataset in {cs.name}")

    def stress(self, cs):
        cs_pairs = cs.from_xyz100(self.xyz_pairs.T).T
        cs_diff = cs_pairs[:, 1] - cs_pairs[:, 0]

        delta = np.sqrt(np.einsum("ij,ij->i", cs_diff, cs_diff))

        alpha = np.dot(self.dist, delta) / np.dot(self.dist, self.dist)
        diff = alpha * self.dist - delta
        return 100 * np.sqrt(np.dot(diff, diff) / np.dot(delta, delta))
