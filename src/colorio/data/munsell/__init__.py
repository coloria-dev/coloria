import json
import pathlib
from typing import Type

import matplotlib.pyplot as plt
import numpy as np

from ...cs import XYY, ColorSpace
from ...illuminants import whitepoints_cie1931
from ..helpers import create_cs_class_instance

this_dir = pathlib.Path(__file__).resolve().parent


class Munsell:
    def __init__(self):
        with open(this_dir / "real.json") as f:
            data = json.load(f)

        self.h = np.array(data["h"])
        self.V = np.array(data["V"])
        self.C = np.array(data["C"])
        xyy100 = np.array([data["x"], data["y"], data["Y"]])
        self.xyz100 = XYY(100).to_xyz100(xyy100)

        # Whitepoint and CIECAM02 info from the JzAzBz paper:
        self.whitepoint_xyz100 = whitepoints_cie1931["C"]
        self.L_A = 64
        self.c = 0.69
        self.Y_b = 20

        with open(this_dir / "lightness.json") as f:
            self.lightness = json.load(f)

    def plot(self, cs: ColorSpace, V: int):
        # pick the data from the given munsell level
        xyz100 = self.xyz100[:, V == self.V]
        pts = cs.from_xyz100(xyz100)

        rgb = cs.to_rgb1(pts)

        # plot the ones that cannot be represented in SRGB in black
        is_legal_srgb = np.all((0 <= rgb) & (rgb <= 1), axis=0)
        fill = rgb.T.copy()
        fill[~is_legal_srgb] = [1.0, 1.0, 1.0]
        edge = rgb.T.copy()
        edge[~is_legal_srgb] = [0.0, 0.0, 0.0]

        idx = [0, 1, 2]
        assert cs.k0 is not None
        k1, k2 = idx[: cs.k0] + idx[cs.k0 + 1 :]

        plt.scatter(pts[k1], pts[k2], marker="o", color=fill, edgecolors=edge)

        plt.title(f"Munsell points at lightness V={V} in {cs.name}")
        plt.xlabel(cs.labels[k1])
        plt.ylabel(cs.labels[k2], rotation=0)
        plt.axis("equal")
        return plt

    def plot_lightness(self, cs: ColorSpace):
        L0_ = cs.from_xyz100(np.zeros(3))[cs.k0]
        L_ = cs.from_xyz100(self.xyz100)[cs.k0] - L0_
        ref = self.V
        alpha = np.dot(ref, L_) / np.dot(ref, ref)

        # plot lightness curve
        v, y = self.lightness
        v = np.asarray(v)
        plt.plot(y, alpha * v, label="scaled Munsell lightness")
        # plt.grid()
        plt.xlabel("Y")
        plt.ylabel(cs.labels[cs.k0], rotation=0)
        plt.title(f"{cs.name} lightness of Munsell samples")

        y_vals = []
        l_avg2 = []
        l_err0 = []
        l_err1 = []

        y_vals2 = []
        l_vals = []
        for k in range(1, 10):
            idx = self.V == k
            y_vals.append(self.xyz100[1, idx][0])
            avg2 = np.sqrt(np.mean((L_[idx] ** 2)))
            # avg2 = np.mean(L_[idx])
            l_avg2.append(avg2)
            l_err0.append(avg2 - np.min(L_[idx]))
            l_err1.append(np.max(L_[idx]) - avg2)
            #
            y_vals2.append(self.xyz100[1, idx])
            l_vals.append(L_[idx])

        plt.errorbar(
            y_vals,
            l_avg2,
            yerr=[l_err0, l_err1],
            fmt="o",
            label=f"{cs.name} lightness",
        )

        # plt.scatter(
        #     np.concatenate(y_vals2),
        #     np.concatenate(l_vals),
        #     color="C1",
        #     marker="o",
        #     label=f"{cs.name} lightness",
        # )
        return plt

    def stress_lightness(self, cs_class: Type[ColorSpace]) -> float:
        cs = create_cs_class_instance(
            cs_class, self.whitepoint_xyz100, self.c, self.Y_b, self.L_A
        )

        ref = self.V

        # Move L0 into origin for translation invariance
        assert cs.k0 is not None
        L0_ = cs.from_xyz100(np.zeros(3))[cs.k0]
        L_ = cs.from_xyz100(self.xyz100)[cs.k0] - L0_

        alpha = np.dot(ref, L_) / np.dot(ref, ref)
        diff = alpha * ref - L_
        return 100 * np.sqrt(np.dot(diff, diff) / np.dot(L_, L_))

    stress = stress_lightness
