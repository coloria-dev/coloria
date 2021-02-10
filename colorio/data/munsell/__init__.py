import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np

from ...cs import XYY

this_dir = pathlib.Path(__file__).resolve().parent


class Munsell:
    def __init__(self):
        with open(this_dir / "real.json") as f:
            data = json.load(f)

        self.h = np.array(data["h"])
        self.V = np.array(data["V"])
        self.C = np.array(data["C"])
        self.xyy100 = np.array([data["x"], data["y"], data["Y"]])

        with open(this_dir / "lightness.json") as f:
            self.lightness = json.load(f)

    def show(self, *args, **kwargs):
        self.plot(*args, **kwargs)
        plt.show()

    def savefig(self, filename, *args, **kwargs):
        plt.figure()
        self.plot(*args, **kwargs)
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        plt.close()

    def plot(self, cs, V):
        # pick the data from the given munsell level
        xyy = self.xyy100[:, V == self.V]

        x, y, Y = xyy
        xyz100 = np.array([Y / y * x, Y, Y / y * (1 - x - y)])
        pts = cs.from_xyz100(xyz100)

        rgb = cs.to_rgb1(pts)
        is_legal_srgb = np.all((0 <= rgb) & (rgb <= 1), axis=0)

        # plot the ones that cannot be represented in SRGB in black
        is_legal_srgb = np.all((0 <= rgb) & (rgb <= 1), axis=0)
        fill = rgb.T.copy()
        fill[~is_legal_srgb] = [1.0, 1.0, 1.0]
        edge = rgb.T.copy()
        edge[~is_legal_srgb] = [0.0, 0.0, 0.0]

        idx = [0, 1, 2]
        k1, k2 = idx[: cs.k0] + idx[cs.k0 + 1 :]

        plt.scatter(pts[k1], pts[k2], marker="o", color=fill, edgecolors=edge)

        plt.title(f"Munsell points at lightness V={V} in {cs.name}")
        plt.xlabel(cs.labels[k1])
        plt.ylabel(cs.labels[k2], rotation=0)
        plt.axis("equal")

    def show_lightness(self, *args, **kwargs):
        self.plot_lightness(*args, **kwargs)
        plt.show()

    def savefig_lightness(self, filename, *args, **kwargs):
        self.plot_lightness(*args, **kwargs)
        plt.savefig(filename, transparent=True, bbox_inches="tight")

    def plot_lightness(self, cs):
        # print(self.xyy100.T)
        # exit(1)

        L0_ = cs.from_xyz100(np.zeros(3))[cs.k0]
        L_ = cs.from_xyz100(XYY(100).to_xyz100(self.xyy100))[cs.k0] - L0_
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
        l_avg = []
        l_err0 = []
        l_err1 = []

        y_vals2 = []
        l_vals = []
        for k in range(1, 10):
            idx = self.V == k
            y_vals.append(self.xyy100[2, idx][0])
            avg = np.average(L_[idx])
            l_avg.append(avg)
            l_err0.append(avg - np.min(L_[idx]))
            l_err1.append(np.max(L_[idx]) - avg)
            #
            y_vals2.append(self.xyy100[2, idx])
            l_vals.append(L_[idx])

        plt.errorbar(
            y_vals,
            l_avg,
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

    def stress_lightness(self, cs):
        ref = self.V

        # Move L0 into origin for translation invariance
        L0_ = cs.from_xyz100(np.zeros(3))[cs.k0]
        L_ = cs.from_xyz100(XYY(100).to_xyz100(self.xyy100))[cs.k0] - L0_

        alpha = np.dot(ref, L_) / np.dot(ref, ref)
        diff = alpha * ref - L_
        return 100 * np.sqrt(np.dot(diff, diff) / np.dot(L_, L_))
