"""
David L. MacAdam,
Uniform color scales,
Journal of the Optical Society of America, Vol. 64, Issue 12, pp. 1691-1702,
1974,
<https://doi.org/10.1364/JOSA.64.001691>.
"""
import json
import pathlib
from typing import Type

import matplotlib.pyplot as plt
import numpy as np

from ...cs import ColorCoordinates, ColorSpace, convert
from ...illuminants import whitepoints_cie1964
from ..color_distance import ColorDistanceDataset
from ..helpers import create_cs_class_instance


class MacAdam1974(ColorDistanceDataset):
    def __init__(self):
        # Extract ellipse centers and offsets from MacAdams data
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "table2.json") as f:
            data = json.load(f)

        # From the article:
        #
        # > Instead, the committee decided to study a limited sample, consisting of 43
        # > colors (made in the form of 5-cm-hexagonal, matte-finish, painted ceramic
        # > tiles, all having approximately the same luminous reflectance (30%) in CIE
        # > Standard D65 daylight, [...]
        #
        # > The CIE specifications for D65 and for the 1964 supplementary observer for
        # > 100 visual field, based on the agreed-upon spectral reflectance data, are
        # > given in Table II.
        #
        # The angle of a 5cm tile viewed with 50 cm distance is
        #
        #   arctan(5cm / 40cm) ~= 7.1 degrees
        #
        # Close enough to the 10-degree-observer I guess?
        #
        self.whitepoint_xyz100 = whitepoints_cie1964["D65"]
        # Assume the surround parameters of a light booth
        self.c = 0.69
        self.Y_b = 20
        self.L_A = 60

        t = dict(zip(data.keys(), range(len(data))))
        xyy100_tiles = ColorCoordinates(
            np.array([[val[0], val[1], val[2]] for val in data.values()]).T, "xyy100"
        )
        self.xyz100_tiles = convert(xyy100_tiles, "xyz100")

        with open(this_dir / "table1.json") as f:
            data = json.load(f)

        d = np.array([item[3] for item in data])
        pairs = np.array([[t[item[1]], t[item[2]]] for item in data])

        # for plotting, only consider the first 43 tiles which are all of approximately
        # the same lightness
        # TODO plot 3D tetrahedral pairs, too
        self.is_flat_pair = np.all(pairs <= 43, axis=1)

        super().__init__("MacAdam (1974)", d, self.xyz100_tiles.data.T[pairs])

    def plot(self, cs_class: Type[ColorSpace]):
        cs = create_cs_class_instance(
            cs_class, self.whitepoint_xyz100, self.c, self.Y_b, self.L_A
        )

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
        d = self.target_dist[self.is_flat_pair]
        alpha = np.dot(d, delta) / np.dot(d, d)
        d *= alpha

        keep = [True, True, True]
        keep[cs.k0] = False

        # make the arrows and dots the same color as the axes labels
        color = plt.gca().xaxis.label.get_color()

        # plot arrows
        for target_dist, pair in zip(d, pairs):
            # arrow from pair[0] to pair[1]
            base = pair[0][keep]
            diff = pair[1] - pair[0]
            v = diff / np.linalg.norm(diff, 2) * target_dist / 2
            v = v[keep]
            plt.plot(
                [base[0], base[0] + v[0]],
                [base[1], base[1] + v[1]],
                color=color,
                linewidth=3.0,
                alpha=0.3,
            )
            # arrow from pair[1] to pair[0]
            base = pair[1][keep]
            v = -v
            plt.plot(
                [base[0], base[0] + v[0]],
                [base[1], base[1] + v[1]],
                color=color,
                linewidth=3.0,
                alpha=0.3,
            )

        # plot colors dots for the first 43 tiles
        tiles = ColorCoordinates(self.xyz100_tiles.data[:, :43], "xyz100")
        coords = convert(tiles, cs)

        plt.scatter(
            coords.hue[0],
            coords.hue[1],
            marker="s",
            color=convert(tiles, "srgb1", mode="clip").data.T,
            edgecolors="w",
            zorder=2,
        )

        plt.gca().set_aspect("equal")
        plt.title(f"MacAdam 1974 color distance data for {cs.name}")

        labels = cs.hue_labels
        plt.xlabel(labels[0])
        plt.ylabel(labels[1], rotation=0)
        return plt
