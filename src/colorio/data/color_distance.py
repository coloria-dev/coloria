from __future__ import annotations

from typing import Callable, Type

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike

from ..cs import CIELAB, ColorSpace
from .helpers import create_cs_class_instance, stress_absolute, stress_relative


class ColorDistanceDataset:
    def __init__(
        self,
        name: str,
        target_dist: ArrayLike,
        xyz_pairs: ArrayLike,
        weights: float | ArrayLike = 1.0,
    ):
        self.name = name
        self.target_dist = np.asarray(target_dist)
        n = len(self.target_dist)
        self.xyz_pairs = np.asarray(xyz_pairs)
        assert n == len(self.xyz_pairs)
        self.weights = np.asarray(weights)

    def plot(self, cs_class: Type[ColorSpace]):
        cs = create_cs_class_instance(
            cs_class, self.whitepoint_xyz100, self.c, self.Y_b, self.L_A
        )

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

    # TODO Python 3.8+: variant: Literal["absolute"] | Literal["relative"]
    def stress(self, cs_class: Type[ColorSpace], variant: str = "absolute") -> float:
        cs = create_cs_class_instance(
            cs_class, self.whitepoint_xyz100, self.c, self.Y_b, self.L_A
        )

        # compute Euclidean distance in the given colorspace
        cs_pairs = cs.from_xyz100(self.xyz_pairs.T).T
        cs_diff = cs_pairs[:, 1] - cs_pairs[:, 0]
        delta = np.sqrt(np.einsum("ij,ij->i", cs_diff, cs_diff))
        fun = stress_absolute if variant == "absolute" else stress_relative
        return fun(self.target_dist, delta, self.weights)

    def stress_lab_diff(self, fun: Callable, variant: str = "absolute") -> float:
        """Same as stress(), but you can provide a color difference function that
        receives two LAB values and returns their scalar distance.
        """
        cielab = CIELAB(self.whitepoint_xyz100)
        lab_pairs = cielab.from_xyz100(self.xyz_pairs.T)
        delta = fun(lab_pairs[:, 0], lab_pairs[:, 1])
        fun = stress_absolute if variant == "absolute" else stress_relative
        return fun(self.target_dist, delta, self.weights)
