"""
David L. MacAdam,
Visual Sensitivities to Color Differences in Daylight,
Journal of the Optional Society of America,
Volume 32, May, 1942, Number 5,
https://doi.org/10.1364/JOSA.32.000247
"""
import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np

from ...cs import XYY
from ..helpers import Dataset


class MacAdam1942(Dataset):
    def __init__(self):
        # Extract ellipse centers and offsets from MacAdams data
        this_dir = pathlib.Path(__file__).resolve().parent
        with open(this_dir / "table3.json") as f:
            data = json.load(f)
        #
        # collect the ellipse centers and offsets
        self.xy_centers = []
        self.xy_offsets = []
        for datak in data:
            # collect ellipse points
            _, _, _, _, delta_y_delta_x, delta_s = np.array(datak["data"]).T
            offset = (
                np.array([np.ones(delta_y_delta_x.shape[0]), delta_y_delta_x])
                / np.sqrt(1 + delta_y_delta_x ** 2)
                * delta_s
            )
            if offset.shape[1] < 2:
                continue
            self.xy_centers.append(np.array([datak["x"], datak["y"]]))
            self.xy_offsets.append(np.column_stack([+offset, -offset]))

    def plot(self, cs, ellipse_scaling=10.0):
        Y = 50.0
        xyy100_centers = []
        xyy100_points = []
        for c, off in zip(self.xy_centers, self.xy_offsets):
            xyy100_centers.append(np.array([*c, Y]))
            p = (c + off.T).T
            xyy100_points.append(np.array([*p, np.full(p.shape[1], Y)]))

        _plot_ellipses(xyy100_centers, xyy100_points, cs, ellipse_scaling)
        plt.title(f"MacAdam ellipses for {cs.name}")

        # cs.plot_visible_slice(
        #     lightness,
        #     outline_prec=outline_prec,
        #     fill_color=visible_gamut_fill_color,
        # )
        # if plot_srgb_gamut:
        #     cs.plot_rgb_slice(lightness)

    def stress(self, cs, Y100: float) -> float:
        xyy100 = XYY(100)

        dists = []
        for c, off in zip(self.xy_centers, self.xy_offsets):
            pts = (c + off.T).T

            # get ellipse center in transformed space
            xyY1_center = np.append(c, Y100)
            c = cs.from_xyz100(xyy100.to_xyz100(xyY1_center))

            # get ellipse points in transformed space
            xyY1_pts = np.array([*pts, np.full(pts.shape[1], Y100)])
            pts = cs.from_xyz100(xyy100.to_xyz100(xyY1_pts))

            # compute the distance in the transformed space
            diff = (c - pts.T).T
            dists.append(np.sqrt(np.einsum("ij,ij->j", diff, diff)))

        delta = np.concatenate(dists)

        diff = np.average(delta) - delta
        return 100 * np.sqrt(np.dot(diff, diff) / np.dot(delta, delta))
