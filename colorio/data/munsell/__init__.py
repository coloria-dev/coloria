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

    def stress_lightness(self, cs):
        # # Each level (1-9) in Munsell is associated with a fixed Y value. Associated
        # # them with another.
        # d = {}
        # d[0] = 0.0
        # for k in range(1, 10):
        #     idx = np.searchsorted(self.V, k)
        #     d[k] = self.xyy100[2, idx]
        # d[10] = 102.5
        # # TODO 406 S. M. NEWHALL, D. NICKERSON, AND D. B. JUDD

        # print(d)
        # exit(1)

        # Move L0 into origin for translation invariance
        L0_ = cs.from_xyz100(np.zeros(3))[cs.k0]
        L_ = cs.from_xyz100(XYY(100).to_xyz100(self.xyy100))[cs.k0] - L0_

        alpha = np.dot(self.V, L_) / np.dot(self.V, self.V)
        diff = alpha * self.V - L_
        return 100 * np.sqrt(np.dot(diff, diff) / np.dot(L_, L_))
