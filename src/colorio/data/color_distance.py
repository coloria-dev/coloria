from typing import Callable, Optional

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike

from ..cs import CIELAB, ColorSpace


class ColorDistanceDataset:
    def __init__(
        self, name: str, dist, xyz_pairs: ArrayLike, weights: Optional[ArrayLike] = None
    ):
        self.name = name
        n = len(dist)
        self.xyz_pairs = np.asarray(xyz_pairs)
        assert n == len(self.xyz_pairs)
        self.dist = np.asarray(dist)
        self.weights = np.ones(n) if weights is None else np.asarray(weights)

    def plot(self, cs: ColorSpace):
        coords = cs.from_xyz100(self.xyz_pairs.T).T

        # reorder the coords such that the lightness in the last (the z-)component
        assert cs.k0 is not None
        coords = np.roll(coords, 2 - cs.k0, axis=0)
        labels = np.roll(cs.labels, 2 - cs.k0, axis=0)

        ax = plt.axes(projection="3d")
        for pair in coords:
            ax.plot(pair[:, 0], pair[:, 1], pair[:, 2], "-k")

        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_zlabel(labels[2])
        ax.set_title(f"{self.name} dataset in {cs.name}")
        return plt

    def stress(self, cs: ColorSpace, variant: str = "absolute"):
        # compute Euclidean distance in colorspace cs
        cs_pairs = cs.from_xyz100(self.xyz_pairs.T).T
        cs_diff = cs_pairs[:, 1] - cs_pairs[:, 0]
        delta = np.sqrt(np.einsum("ij,ij->i", cs_diff, cs_diff))
        return self._stress(delta, variant)

    def stress_lab_diff(self, fun: Callable, variant: str = "absolute"):
        """Same as stress(), but you can provide a color difference function that
        receives two LAB values and returns their scalar distance.
        """
        lab_pairs = CIELAB().from_xyz100(self.xyz_pairs.T)
        delta = fun(lab_pairs[:, 0], lab_pairs[:, 1])
        return self._stress(delta, variant)

    def _stress(self, delta: ArrayLike, variant: str):
        delta = np.asarray(delta)
        if variant == "absolute":
            # regular old stress
            wdist = self.weights * self.dist
            alpha = np.dot(wdist, delta) / np.dot(wdist, self.dist)
            diff = alpha * self.dist - delta
            val = np.dot(self.weights * diff, diff) / np.dot(
                self.weights * delta, delta
            )
        else:
            assert variant == "relative", f"Illegal variant {variant}."
            alpha = np.sum(self.weights * self.dist) / np.sum(
                self.weights * self.dist ** 2 / delta
            )
            diff = alpha * self.dist - delta
            val = np.sum(self.weights * diff ** 2 / delta) / np.sum(
                self.weights * delta
            )
        return 100 * np.sqrt(val)
