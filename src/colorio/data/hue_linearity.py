import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike

from ..cs import ColorSpace


class HueLinearityDataset:
    def __init__(self, name: str, whitepoint_xyz100: ArrayLike, arms):
        self.name = name
        self.whitepoint_xyz100 = np.asarray(whitepoint_xyz100)
        self.arms = arms

    def plot(self, cs: ColorSpace):
        # k0 is the coordinate that corresponds to "lightness"
        assert cs.k0 is not None
        no_lightness = [True, True, True]
        no_lightness[cs.k0] = False

        wp = cs.from_xyz100(self.whitepoint_xyz100)[no_lightness]
        all_pts = []
        all_rgb1 = []
        for xyz in self.arms:
            pts = cs.from_xyz100(xyz)
            rgb1 = cs.to_rgb1(pts)
            pts = pts[no_lightness]

            # get the eigenvector corresponding to the larger eigenvalue
            pts_wp = (pts.T - wp).T
            vals, vecs = np.linalg.eigh(pts_wp @ pts_wp.T)
            v = vecs[:, 0] if vals[0] > vals[1] else vecs[:, 1]
            # invert if necessary
            if np.dot(v, np.average(pts_wp, axis=1)) < 0:
                v = -v

            length = np.max(np.linalg.norm(pts_wp, axis=0))
            end_point = wp + length * v
            plt.plot(
                [wp[0], end_point[0]], [wp[1], end_point[1]], "-", color="0.5", zorder=0
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

        l0, l1 = cs.labels[no_lightness]
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

    def stress(self, cs: ColorSpace):
        """Compute the TLS residuals for each of the arms."""
        # return _compute_straight_line_stress(cs, self.whitepoint, self.arms)
        # def _compute_straight_line_stress(cs, wp, d):
        # remove the row corresponding to lightness
        assert cs.k0 is not None
        idx = [True, True, True]
        idx[cs.k0] = False

        wp_cs = cs.from_xyz100(self.whitepoint_xyz100)[idx]

        s2 = []
        for dd in self.arms:
            vals = cs.from_xyz100(dd)[idx]
            # move values such that whitepoint is in the origin
            vals = (vals.T - wp_cs).T
            # could also be computed explicitly
            s_max, s_min = np.linalg.svd(vals, compute_uv=False)
            s2.append(s_min / s_max)
            # plt.plot(vals[0], vals[1], "x")
            # plt.gca().set_aspect("equal")
            # plt.show()
        return 100 * np.array(s2)
