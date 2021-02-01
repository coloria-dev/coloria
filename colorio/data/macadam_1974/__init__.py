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
from ..helpers import ColorDistanceDataset


class MacAdam1974(ColorDistanceDataset):
    def __init__(self):
        # Extract ellipse centers and offsets from MacAdams data
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "table2.json") as f:
            data = json.load(f)

        t = dict(zip(data.keys(), range(len(data))))
        xyy1_tiles = np.array([[val[0], val[1], val[2]] for val in data.values()])
        xyz100_tiles = XYY(100).to_xyz100(xyy1_tiles.T).T

        with open(this_dir / "table1.json") as f:
            data = json.load(f)

        d = np.array([item[3] for item in data])
        pairs = np.array([[t[item[1]], t[item[2]]] for item in data])

        # for plotting, only consider the first 43 tiles which are all of approximately
        # the same lightness
        # TODO plot 3D tetrahedral pairs, too
        self.is_flat_pair = np.all(pairs <= 43, axis=1)

        super().__init__("MacAdam (1974)", d, xyz100_tiles[pairs])

    def plot(self, cs):
        pairs = self.xyz_pairs[self.is_flat_pair]
        pairs = cs.from_xyz100(pairs.T).T

        # print(diff.shape)
        # exit(1)
        # d = self.d[self.is_flat_pair]
        # pairs = self.pairs[self.is_flat_pair]
        # xyz100_tiles = self.xyz100_tiles[:, :43]
        # pts = cs.from_xyz100(xyz100_tiles)

        # # Plot the tile points.
        # pts_2d = np.delete(pts, cs.k0, axis=0)
        # plt.plot(pts_2d[0], pts_2d[1], "ok", fillstyle="none")

        # for k, pt in enumerate(pts.T):
        #     # plt.text(pt[0], pt[1], k + 1)
        #     plt.plot(pt[0], pt[1], "ok", fillstyle="none")

        # scale the distances
        delta = np.linalg.norm(pairs[:, 0] - pairs[:, 1], axis=1)
        d = self.dist[self.is_flat_pair]
        alpha = np.dot(d, delta) / np.dot(d, d)
        d *= alpha

        keep = [True, True, True]
        keep[cs.k0] = False
        # plot arrows
        for target_dist, pair in zip(d, pairs):
            # arrow from pair[0] to pair[1]
            base = pair[0][keep]
            diff = pair[1] - pair[0]
            v = diff / np.linalg.norm(diff, 2) * target_dist / 2
            v = v[keep]
            plt.arrow(
                base[0], base[1], v[0], v[1], length_includes_head=True, color="k"
            )
            # arrow from pair[1] to pair[0]
            base = pair[1][keep]
            v = -v
            plt.arrow(
                base[0], base[1], v[0], v[1], length_includes_head=True, color="k"
            )

        # remove lightness coord
        pairs = pairs[..., keep]
        plt.scatter(pairs[..., 0], pairs[..., 1], color="k")

        plt.gca().set_aspect("equal")
        plt.title(f"MacAdam 1974 color distance data for {cs.name}")

        labels = cs.labels[keep]
        plt.xlabel(labels[0])
        plt.ylabel(labels[1], rotation=0)
