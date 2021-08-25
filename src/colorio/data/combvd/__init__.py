from typing import Callable, Type

import numpy as np

from ...cs import CIELAB, ColorSpace
from ..bfd_p import BfdP
from ..helpers import create_cs_class_instance, stress_absolute, stress_relative
from ..leeds import Leeds
from ..rit_dupont import RitDupont
from ..witt import Witt


class COMBVD:
    def __init__(self):
        self.bfd_p = BfdP()
        self.leeds = Leeds()
        self.rit_dupont = RitDupont()
        self.witt = Witt()

    def stress(self, cs_class: Type[ColorSpace], variant: str = "absolute"):
        deltas = []
        target_dist = []
        weights = []

        # bfd-p
        for xyz, wp, td in zip(
            self.bfd_p.xyz_pairs, self.bfd_p.whitepoints, self.bfd_p.target_dist
        ):
            cs = create_cs_class_instance(
                cs_class, wp, self.bfd_p.c, self.bfd_p.Y_b, self.bfd_p.L_A
            )
            cs_pairs = cs.from_xyz100(xyz.T).T
            cs_diff = cs_pairs[:, 1] - cs_pairs[:, 0]
            deltas.append(np.sqrt(np.einsum("ij,ij->i", cs_diff, cs_diff)))
            target_dist.append(td)
            weights.append(np.full(len(td), 1.0))

        # leeds, witt, rit_dupont
        for dataset, weight in [
            (self.leeds, 9.0),
            (self.rit_dupont, 9.0),
            (self.witt, 7.0),
        ]:
            cs = create_cs_class_instance(
                cs_class, dataset.whitepoint_xyz100, dataset.c, dataset.Y_b, dataset.L_A
            )
            cs_pairs = cs.from_xyz100(dataset.xyz_pairs.T).T
            cs_diff = cs_pairs[:, 1] - cs_pairs[:, 0]
            deltas.append(np.sqrt(np.einsum("ij,ij->i", cs_diff, cs_diff)))
            target_dist.append(dataset.target_dist)
            weights.append(np.full(len(dataset.target_dist), weight))

        target_dist = np.concatenate(target_dist)
        delta = np.concatenate(deltas)
        weights = np.concatenate(weights)

        fun = stress_absolute if variant == "absolute" else stress_relative
        return fun(target_dist, delta, weights)

    def stress_lab_diff(self, fun: Callable, variant: str = "absolute") -> float:
        deltas = []
        target_dist = []
        weights = []

        # bfd-p
        for xyz, wp, td in zip(
            self.bfd_p.xyz_pairs, self.bfd_p.whitepoints, self.bfd_p.target_dist
        ):
            lab_pairs = CIELAB(wp).from_xyz100(xyz.T)
            deltas.append(fun(lab_pairs[:, 0], lab_pairs[:, 1]))
            target_dist.append(td)
            weights.append(np.full(len(td), 1.0))

        # leeds, witt, rit_dupont
        for dataset, weight in [
            (self.leeds, 9.0),
            (self.rit_dupont, 9.0),
            (self.witt, 7.0),
        ]:
            lab_pairs = CIELAB(dataset.whitepoint_xyz100).from_xyz100(
                dataset.xyz_pairs.T
            )
            deltas.append(fun(lab_pairs[:, 0], lab_pairs[:, 1]))
            target_dist.append(dataset.target_dist)
            weights.append(np.full(len(dataset.target_dist), weight))

        target_dist = np.concatenate(target_dist)
        delta = np.concatenate(deltas)
        weights = np.concatenate(weights)

        fun = stress_absolute if variant == "absolute" else stress_relative
        return fun(target_dist, delta, weights)
