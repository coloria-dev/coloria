"""
David L. MacAdam,
Uniform color scales,
Journal of the Optical Society of America, Vol. 64, Issue 12, pp. 1691-1702,
1974,
<https://doi.org/10.1364/JOSA.64.001691>.
"""
import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np

from ...cs import XYY
from ..helpers import Dataset


class MacAdam1974(Dataset):
    def __init__(self):
        # Extract ellipse centers and offsets from MacAdams data
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "table2.json") as f:
            data = json.load(f)

        t = dict(zip(data.keys(), range(len(data))))
        xyy1_tiles = np.array([[val[0], val[1], val[2]] for val in data.values()])
        self.xyz100_tiles = XYY(100).to_xyz100(xyy1_tiles.T)

        with open(this_dir / "table1.json") as f:
            data = json.load(f)

        self.d = np.array([item[3] for item in data])
        self.pairs = np.array([[t[item[1]], t[item[2]]] for item in data])

    def plot(self, cs):
        # only consider the first 43 tiles which are all of approximately the same lightness
        # TODO plot 3D tetrahedral pairs, too
        d = self.d[np.all(self.pairs <= 43, axis=1)]
        pairs = self.pairs[np.all(self.pairs <= 43, axis=1)]
        xyz100_tiles = self.xyz100_tiles[:, :43]
        pts = cs.from_xyz100(xyz100_tiles)

        # Plot the tile points.
        pts_2d = np.delete(pts, cs.k0, axis=0)
        plt.plot(pts_2d[0], pts_2d[1], "ok", fillstyle="none")

        # for k, pt in enumerate(pts.T):
        #     # plt.text(pt[0], pt[1], k + 1)
        #     plt.plot(pt[0], pt[1], "ok", fillstyle="none")

        # scale the distances
        diff = pts[:, pairs]
        delta = np.linalg.norm(diff[..., 0] - diff[..., 1], axis=0)
        alpha = np.dot(d, delta) / np.dot(d, d)
        d *= alpha

        # plot arrow
        for dist, pair in zip(d, pairs):
            # arrow from pair[0] to pair[1]
            base = pts[:, pair[0]]
            diff = pts[:, pair[1]] - pts[:, pair[0]]
            v = diff / np.linalg.norm(diff, 2) * dist / 2
            base = np.delete(base, cs.k0)
            v = np.delete(v, cs.k0)
            plt.arrow(
                base[0], base[1], v[0], v[1], length_includes_head=True, color="k"
            )
            # arrow from pair[1] to pair[0]
            base = pts[:, pair[1]]
            base = np.delete(base, cs.k0)
            v = -v
            plt.arrow(
                base[0], base[1], v[0], v[1], length_includes_head=True, color="k"
            )

        plt.gca().set_aspect("equal")
        plt.title(f"MacAdam 1974 color distance data for {cs.name}")

    def stress(self, cs):
        pts = cs.from_xyz100(self.xyz100_tiles)

        diff = pts[:, self.pairs]
        delta = np.linalg.norm(diff[..., 0] - diff[..., 1], axis=0)

        alpha = np.dot(self.d, delta) / np.dot(self.d, self.d)
        val = np.dot(alpha * self.d - delta, alpha * self.d - delta) / np.dot(
            delta, delta
        )
        return 100 * np.sqrt(val)
