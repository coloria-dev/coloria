from typing import Callable, Optional

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike

from ..cs import CIELAB, ColorSpace


class ColorDistanceDataset:
    def __init__(
        self, name: str, dist, xyz_pairs: ArrayLike, weights: Optional[ArrayLike] = None
    ):
        self.name = name
        n = len(dist)
        self.xyz_pairs = np.asarray(xyz_pairs)
        assert n == len(self.xyz_pairs)
        self.dist = np.asarray(dist)
        self.weights = np.ones(n) if weights is None else np.asarray(weights)

    def plot(self, cs: ColorSpace):
        coords = cs.from_xyz100(self.xyz_pairs.T).T

        # reorder the coords such that the lightness in the last (the z-)component
        assert cs.k0 is not None
        coords = np.roll(coords, 2 - cs.k0, axis=0)
        labels = np.roll(cs.labels, 2 - cs.k0, axis=0)

        ax = plt.axes(projection="3d")
        for pair in coords:
            ax.plot(pair[:, 0], pair[:, 1], pair[:, 2], "-k")

        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_zlabel(labels[2])
        ax.set_title(f"{self.name} dataset in {cs.name}")
        return plt

    def stress(self, cs: ColorSpace, variant: str = "absolute"):
        # compute Euclidean distance in colorspace cs
        cs_pairs = cs.from_xyz100(self.xyz_pairs.T).T
        cs_diff = cs_pairs[:, 1] - cs_pairs[:, 0]
        delta = np.sqrt(np.einsum("ij,ij->i", cs_diff, cs_diff))
        return self._stress(delta, variant)

    def stress_lab_diff(self, fun: Callable, variant: str = "absolute"):
        """Same as stress(), but you can provide a color difference function that
        receives two LAB values and returns their scalar distance.
        """
        lab_pairs = CIELAB().from_xyz100(self.xyz_pairs.T)
        delta = fun(lab_pairs[:, 0], lab_pairs[:, 1])
        return self._stress(delta, variant)

    def _stress(self, delta: ArrayLike, variant: str):
        delta = np.asarray(delta)
        if variant == "absolute":
            # regular old stress
            wdist = self.weights * self.dist
            alpha = np.dot(wdist, delta) / np.dot(wdist, self.dist)
            diff = alpha * self.dist - delta
            val = np.dot(self.weights * diff, diff) / np.dot(
                self.weights * delta, delta
            )
        else:
            assert variant == "relative", f"Illegal variant {variant}."
            alpha = np.sum(self.weights * self.dist) / np.sum(
                self.weights * self.dist ** 2 / delta
            )
            diff = alpha * self.dist - delta
            val = np.sum(self.weights * diff ** 2 / delta) / np.sum(
                self.weights * delta
            )
        return 100 * np.sqrt(val)


class HueLinearityDataset:
    def __init__(self, name: str, whitepoint_xyz100: ArrayLike, arms):
        self.name = name
        self.whitepoint_xyz100 = np.asarray(whitepoint_xyz100)
        self.arms = arms

    def plot(self, colorspace: ColorSpace):
        # k0 is the coordinate that corresponds to "lightness"
        assert colorspace.k0 is not None
        no_lightness = [True, True, True]
        no_lightness[colorspace.k0] = False

        wp = colorspace.from_xyz100(self.whitepoint_xyz100)[no_lightness]
        all_pts = []
        all_rgb1 = []
        for xyz in self.arms:
            pts = colorspace.from_xyz100(xyz)
            rgb1 = colorspace.to_rgb1(pts)
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

        l0, l1 = colorspace.labels[no_lightness]
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
        plt.title(f"{self.name} hue linearity data for {colorspace.name}")
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


class EllipseDataset:
    def __init__(self, name: str, xyy100_centers, xyy100_points):
        self.name = name
        self.xyy100_centers = xyy100_centers
        self.xyy100_points = xyy100_points

    def stress(self, cs):
        distances = []
        for xyy100_center, xyy100_pts in zip(self.xyy100_centers, self.xyy100_points):
            # append Y
            cs_center = cs.from_xyz100(_xyy100_to_xyz100(xyy100_center))
            cs_ellips = cs.from_xyz100(_xyy100_to_xyz100(xyy100_pts))
            # compute distances to ellipse center
            diff = (cs_center - cs_ellips.T).T
            distances.append(np.sqrt(np.einsum("ij,ij->j", diff, diff)))

        distances = np.concatenate(distances)
        alpha = np.average(distances)
        return 100 * np.sqrt(np.sum((alpha - distances) ** 2) / np.sum(distances ** 2))

    def plot(self, cs, ellipse_scaling=1.0):
        _plot_ellipses(
            cs, self.xyy100_centers, self.xyy100_points, ellipse_scaling=ellipse_scaling
        )
        plt.title(f"{self.name} ellipses for {cs.name}")
        return plt


def _plot_ellipses(cs, xyy100_centers, xyy100_points, ellipse_scaling):
    from matplotlib.patches import Ellipse
    from scipy.optimize import leastsq

    keep = [True, True, True]
    keep[cs.k0] = False

    # make the ellipses the same color as the axes labels
    color = plt.gca().xaxis.label.get_color()

    for center, points in zip(xyy100_centers, xyy100_points):
        # cut off the irrelevant index
        cs_center = cs.from_xyz100(_xyy100_to_xyz100(center))
        cs_points = cs.from_xyz100(_xyy100_to_xyz100(points))

        # project out lightness component
        tcenter = cs_center[keep]
        tvals = cs_points[keep]

        # Given these new transformed vals, find the ellipse that best fits those
        # points
        X = (tvals.T - tcenter).T

        def f_ellipse(a_b_theta):
            a, b, theta = a_b_theta
            sin_t = np.sin(theta)
            cos_t = np.cos(theta)
            return (
                +(a ** 2) * (X[0] * cos_t + X[1] * sin_t) ** 2
                + b ** 2 * (X[0] * sin_t - X[1] * cos_t) ** 2
                - 1.0
            )

        def jac(a_b_theta):
            a, b, theta = a_b_theta
            x0sin = X[0] * np.sin(theta)
            x0cos = X[0] * np.cos(theta)
            x1sin = X[1] * np.sin(theta)
            x1cos = X[1] * np.cos(theta)
            return np.array(
                [
                    +2 * a * (x0cos + x1sin) ** 2,
                    +2 * b * (x0sin - x1cos) ** 2,
                    +(a ** 2) * 2 * (x0cos + x1sin) * (-x0sin + x1cos)
                    + b ** 2 * 2 * (x0sin - x1cos) * (x0cos + x1sin),
                ]
            ).T

        # We need to use some optimization here to find the new ellipses which best
        # fit the modified data.
        (a, b, theta), _ = leastsq(f_ellipse, [1.0, 1.0, 0.0], Dfun=jac)

        # plot the scaled ellipse
        e = Ellipse(
            xy=tcenter,
            width=ellipse_scaling * 2 / a,
            height=ellipse_scaling * 2 / b,
            angle=theta / np.pi * 180,
            # label=label,
        )
        plt.gca().add_patch(e)
        # e.set_alpha(0.5)

        e.set_facecolor(color)

        # plt.plot(*tcenter, "xk")
        # plt.plot(*tvals, "ok")
        # plt.show()

    plt.gca().set_aspect("equal")
    labels = cs.labels[keep]
    plt.xlabel(labels[0])
    plt.ylabel(labels[1], rotation=0)

    # mpl doesn't update axis limits when adding artists,
    # <https://github.com/matplotlib/matplotlib/issues/19290>.
    # Handle it manually for now.
    tcenters = []
    for center, points in zip(xyy100_centers, xyy100_points):
        cs_center = cs.from_xyz100(_xyy100_to_xyz100(center))
        tcenters.append(np.delete(cs_center, cs.k0))
    tcenters = np.asarray(tcenters).T
    xmin = np.min(tcenters[0])
    xmax = np.max(tcenters[0])
    ymin = np.min(tcenters[1])
    ymax = np.max(tcenters[1])
    width = xmax - xmin
    height = ymax - ymin
    plt.xlim(xmin - 0.2 * width, xmax + 0.2 * width)
    plt.ylim(ymin - 0.2 * height, ymax + 0.2 * height)


def _xyy100_to_xyz100(xyy):
    x, y, Y = xyy
    return np.array([Y / y * x, Y, Y / y * (1 - x - y)])
