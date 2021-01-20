import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np

from ...cs import XYY, SrgbLinear

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
        plt.figure()
        self.plot(*args, **kwargs)
        plt.show()
        plt.close()

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
        vals = cs.from_xyz100(xyz100)

        srgb = SrgbLinear()
        rgb = srgb.from_xyz100(xyz100)
        is_legal_srgb = np.all((0 <= rgb) & (rgb <= 1), axis=0)

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

    def residual_lightness(self, cs):
        # Move L0 into origin for translation invariance
        L0_ = cs.from_xyz100(np.zeros(3))[cs.k0]
        L_ = cs.from_xyz100(XYY(100).to_xyz100(self.xyy100))[cs.k0] - L0_

        alpha = np.dot(self.V, L_) / np.dot(self.V, self.V)
        diff = alpha * self.V - L_
        val = np.sqrt(np.dot(diff, diff) / np.dot(L_, L_))
        return val

    def stress_lightness(self, *args):
        return 100 * self.residual_lightness(*args)
