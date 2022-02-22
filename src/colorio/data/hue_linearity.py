from __future__ import annotations

from typing import Type

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike

from ..cs import ColorCoordinates, ColorSpace, convert
from .helpers import create_cs_class_instance


class HueLinearityDataset:
    def __init__(
        self,
        name: str,
        whitepoint_xyz100: ArrayLike,
        arms,
        neutral_gray: ArrayLike | None = None,
    ):
        self.name = name
        self.whitepoint_xyz100 = np.asarray(whitepoint_xyz100)
        self.arms = arms
        self.neutral_gray = (
            self.whitepoint_xyz100 if neutral_gray is None else neutral_gray
        )

    def plot(self, cs_class: Type[ColorSpace]):
        cs = create_cs_class_instance(
            cs_class, self.whitepoint_xyz100, self.c, self.Y_b, self.L_A
        )
        # k0 is the coordinate that corresponds to "lightness"
        ng = convert(ColorCoordinates(self.neutral_gray, "XYZ100"), cs).hue

        all_pts = []
        all_rgb1 = []
        for xyz in self.arms:
            coords = convert(ColorCoordinates(xyz, "XYZ100"), cs)
            rgb1 = convert(coords, "SRGB1", mode="clip").data

            pts = coords.hue

            # get the eigenvector corresponding to the larger eigenvalue
            pts_ng = (pts.T - ng).T
            vals, vecs = np.linalg.eigh(pts_ng @ pts_ng.T)
            v = vecs[:, 0] if vals[0] > vals[1] else vecs[:, 1]
            # invert if necessary
            if np.dot(v, np.average(pts_ng, axis=1)) < 0:
                v = -v

            length = np.max(np.linalg.norm(pts_ng, axis=0))
            end_point = ng + length * v
            plt.plot(
                [ng[0], end_point[0]], [ng[1], end_point[1]], "-", color="0.5", zorder=0
            )

            all_pts.append(pts)
            all_rgb1.append(rgb1)

        # plot all points in one big scatter plot;
        # this helps when exporting the data to pgfplots
        all_pts = np.hstack(all_pts)
        all_rgb1 = np.hstack(all_rgb1)
        is_legal_srgb = np.all((0 <= all_rgb1) & (all_rgb1 <= 1), axis=0)
        fill = all_rgb1.T.copy()
        fill[~is_legal_srgb] = [1.0, 1.0, 1.0]  # white
        edge = all_rgb1.T.copy()
        edge[~is_legal_srgb] = [0.0, 0.0, 0.0]  # black
        plt.scatter(all_pts[0], all_pts[1], marker="o", color=fill, edgecolors=edge)

        l0, l1 = cs.hue_labels
        plt.xlabel(l0)
        plt.ylabel(l1, rotation=0)
        plt.axis("equal")

        # plt.grid()
        plt.grid(False)
        ax = plt.gca()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        plt.title(f"{self.name} hue linearity data for {cs.name}")
        return plt

    def stress(self, cs_class: Type[ColorSpace]) -> np.ndarray:
        """Compute the TLS residuals for each of the arms."""
        cs = create_cs_class_instance(
            cs_class, self.whitepoint_xyz100, self.c, self.Y_b, self.L_A
        )

        # remove the row corresponding to lightness
        assert cs.k0 is not None
        idx = [True, True, True]
        idx[cs.k0] = False

        ng_cs = cs.from_xyz100(self.neutral_gray)[idx]

        s2 = []
        for dd in self.arms:
            vals = cs.from_xyz100(dd)[idx]
            # move values such that whitepoint is in the origin
            vals = (vals.T - ng_cs).T
            # could also be computed explicitly
            s_max, s_min = np.linalg.svd(vals, compute_uv=False)
            s2.append(s_min / s_max)
            # plt.plot(vals[0], vals[1], "x")
            # plt.gca().set_aspect("equal")
            # plt.show()
        return 100 * np.array(s2)
