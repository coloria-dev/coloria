from typing import Optional, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from .._helpers import _find_Y
from .._tools import get_mono_outline_xy, spectrum_to_xyz100
from ..observers import cie_1931_2
from ._hdr import HdrLinear
from ._srgb import SrgbLinear


class ColorSpace:
    def __init__(
        self,
        name: str,
        labels: Tuple[str, str, str],
        k0: Optional[int],
        is_origin_well_defined: bool = True,
    ):
        self.name = name
        self.labels = np.asarray(labels)
        self.k0 = k0  # the index that corresponds to luminosity
        self.is_origin_well_defined = is_origin_well_defined
        self.srgb_linear = SrgbLinear()

    def __repr__(self):
        return f"<colorio color space {self.name}>"

    def to_xyz100(self, cs_coords):
        raise NotImplementedError("ColorSpace needs to implement to_xyz100()")

    def from_xyz100(self, cs_coords):
        raise NotImplementedError("ColorSpace needs to implement from_xyz100()")

    def to_rgb_linear(self, cs_coords):
        return self.srgb_linear.from_xyz100(self.to_xyz100(cs_coords))

    def to_rgb1(self, cs_coords):
        return self.srgb_linear.to_rgb1(self.to_rgb_linear(cs_coords))

    def to_rgb255(self, cs_coords):
        return self.srgb_linear.to_rgb255(self.to_rgb_linear(cs_coords))

    def from_rgb_linear(self, rgb_lin):
        return self.from_xyz100(self.srgb_linear.to_xyz100(rgb_lin))

    def from_rgb1(self, rgb1):
        return self.from_rgb_linear(self.srgb_linear.from_rgb1(rgb1))

    def from_rgb255(self, rgb255):
        return self.from_rgb_linear(self.srgb_linear.from_rgb255(rgb255))

    def save_visible_gamut(self, observer, illuminant, filename):
        import meshio
        from scipy.spatial import ConvexHull

        lmbda, illu = illuminant
        values = []

        # Iterate over every possible illuminant input and store it in values
        n = len(lmbda)
        values = np.empty((n * (n - 1) + 2, 3))
        k = 0

        # No light
        data = np.zeros(n)
        values[k] = spectrum_to_xyz100((lmbda, illu * data), observer=observer)
        k += 1
        # frequency blocks
        for width in range(1, n):
            data = np.zeros(n)
            data[:width] = 1.0
            for _ in range(n):
                values[k] = spectrum_to_xyz100((lmbda, illu * data), observer=observer)
                k += 1
                data = np.roll(data, shift=1)
        # Full illuminant
        data = np.ones(len(lmbda))
        values[k] = spectrum_to_xyz100((lmbda, illu * data), observer=observer)
        k += 1

        # scale the values such that the Y-coordinate of the white point (last entry)
        # has value 100.
        values *= 100 / values[-1][1]

        cells = ConvexHull(values).simplices

        if not self.is_origin_well_defined:
            values = values[1:]
            cells = cells[~np.any(cells == 0, axis=1)]
            cells -= 1

        pts = self.from_xyz100(values.T).T

        meshio.write_points_cells(filename, pts, cells={"triangle": cells})

    def save_rgb_gamut(self, filename: str, variant: str = "srgb", n: int = 50):
        import meshio
        import meshzoo

        if variant.lower() in ["srgb", "rec709"]:
            rgb_linear = SrgbLinear()
        else:
            assert variant.lower() in ["hdr", "rec2020", "rec2100"]
            rgb_linear = HdrLinear()

        points, cells = meshzoo.cube(nx=n, ny=n, nz=n)

        if not self.is_origin_well_defined:
            # cut off [0, 0, 0] to avoid division by 0 in the xyz conversion
            points = points[1:]
            cells = cells[~np.any(cells == 0, axis=1)]
            cells -= 1

        pts = self.from_rgb_linear(points.T).T
        # pts = self.from_xyz100(rgb_linear.to_xyz100(points.T)).T
        assert pts.shape[1] == 3
        rgb = rgb_linear.to_rgb1(points)
        meshio.write_points_cells(
            filename, pts, {"tetra": cells}, point_data={"srgb": rgb}
        )

    def save_cone_gamut(self, filename, observer, max_Y):
        import meshio
        import pygmsh

        with pygmsh.geo.Geometry() as geom:
            max_stepsize = 4.0e-2
            xy, _ = get_mono_outline_xy(observer, max_stepsize=max_stepsize)

            # append third component
            xy = np.column_stack([xy, np.full(xy.shape[0], 1.0e-5)])

            # Draw a cross.
            poly = geom.add_polygon(xy, mesh_size=max_stepsize)

            axis = [0, 0, max_Y]

            geom.extrude(poly, translation_axis=axis)

            mesh = geom.generate_mesh(verbose=False)
        # meshio.write(filename, mesh)

        pts = self.from_xyz100(_xyy_to_xyz100(mesh.points.T)).T
        meshio.write_points_cells(
            filename, pts, {"tetra": mesh.get_cells_type("tetra")}
        )

    def show_visible_slice(self, *args, **kwargs):
        plt.figure()
        self.plot_visible_slice(*args, **kwargs)
        plt.show()
        plt.close()

    def save_visible_slice(self, filename, *args, **kwargs):
        plt.figure()
        self.plot_visible_slice(*args, **kwargs)
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        plt.close()

    def plot_visible_slice(self, lightness, outline_prec=1.0e-2, fill_color="0.8"):
        # first plot the monochromatic outline
        mono_xy, conn_xy = get_mono_outline_xy(
            observer=cie_1931_2(), max_stepsize=outline_prec
        )

        mono_vals = np.array([_find_Y(self, xy, lightness) for xy in mono_xy])
        conn_vals = np.array([_find_Y(self, xy, lightness) for xy in conn_xy])

        k1, k2 = [k for k in [0, 1, 2] if k != self.k0]
        plt.plot(mono_vals[:, k1], mono_vals[:, k2], "-", color="k")
        plt.plot(conn_vals[:, k1], conn_vals[:, k2], ":", color="k")
        #
        if fill_color is not None:
            xyz = np.vstack([mono_vals, conn_vals[1:]])
            plt.fill(xyz[:, k1], xyz[:, k2], facecolor=fill_color, zorder=0)

        plt.axis("equal")
        plt.xlabel(self.labels[k1])
        plt.ylabel(self.labels[k2])
        plt.title(f"{self.labels[self.k0]} = {lightness}")

    def show_rgb_slice(self, *args, **kwargs):
        plt.figure()
        self.plot_rgb_slice(*args, **kwargs)
        plt.show()
        plt.close()

    def save_rgb_slice(self, filename, *args, **kwargs):
        plt.figure()
        self.plot_rgb_slice(*args, **kwargs)
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        plt.close()

    def plot_rgb_slice(self, lightness: float, n: int = 50, variant: str = "srgb"):
        import meshzoo

        # TODO HDR
        assert variant in ["srgb", "rec709"]

        # Get all RGB values that sum up to 1.
        srgb_vals, triangles = meshzoo.triangle(n=n)
        srgb_vals = srgb_vals.T

        # Use bisection to
        srgb_linear = SrgbLinear()
        tol = 1.0e-5
        # Use zeros() instead of empty() here to avoid invalid values when setting up
        # the cmap below.
        self_vals = np.zeros((srgb_vals.shape[0], 3))
        srgb_linear_vals = np.zeros((srgb_vals.shape[0], 3))
        mask = np.ones(srgb_vals.shape[0], dtype=bool)
        for k, val in enumerate(srgb_vals):
            alpha_min = 0.0
            xyz100 = srgb_linear.to_xyz100(val * alpha_min)
            self_val_min = self.from_xyz100(xyz100)[self.k0]
            if self_val_min > lightness:
                mask[k] = False
                continue

            alpha_max = 1.0 / np.max(val)

            xyz100 = srgb_linear.to_xyz100(val * alpha_max)
            self_val_max = self.from_xyz100(xyz100)[self.k0]
            if self_val_max < lightness:
                mask[k] = False
                continue

            # bisection
            while True:
                alpha = (alpha_max + alpha_min) / 2
                srgb_linear_vals[k] = val * alpha
                xyz100 = srgb_linear.to_xyz100(srgb_linear_vals[k])
                self_val = self.from_xyz100(xyz100)
                if abs(self_val[self.k0] - lightness) < tol:
                    break
                elif self_val[self.k0] < lightness:
                    alpha_min = alpha
                else:
                    assert self_val[self.k0] > lightness
                    alpha_max = alpha
            self_vals[k] = self_val

        # Remove all triangles which have masked corner points
        tri_mask = np.all(mask[triangles], axis=1)
        if ~np.any(tri_mask):
            return
        triangles = triangles[tri_mask]

        # Unfortunately, one cannot use tripcolors with explicit RGB specification (see
        # <https://github.com/matplotlib/matplotlib/issues/10265>). As a workaround,
        # associate range(n) data with the points and create a colormap that associates
        # the integer values with the respective RGBs.
        z = np.arange(srgb_vals.shape[0])
        rgb = srgb_linear.to_rgb1(srgb_linear_vals)
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "gamut", rgb, N=len(rgb)
        )

        k1, k2 = [k for k in [0, 1, 2] if k != self.k0]

        plt.tripcolor(
            self_vals[:, k1],
            self_vals[:, k2],
            triangles,
            z,
            shading="gouraud",
            cmap=cmap,
        )
        # plt.triplot(self_vals[:, k1], self_vals[:, k2], triangles=triangles)

    def show_srgb1_gradient(self, *args, **kwargs):
        plt.figure()
        self.plot_srgb1_gradient(*args, **kwargs)
        plt.show()
        plt.close()

    def plot_srgb1_gradient(self, srgb0, srgb1, n=256):
        srgb = self.get_srgb1_gradient(srgb0, srgb1, n=n)

        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("empty", srgb, n)

        gradient = np.linspace(0.0, 1.0, n)
        gradient = np.vstack((gradient, gradient))
        plt.imshow(gradient, aspect="auto", cmap=cmap)
        plt.axis("off")
        plt.title(f"SRGB gradient in {self.name}")

    def get_srgb1_gradient(self, srgb0, srgb1, n):
        # convert to colorspace
        cs = [self.from_rgb1(srgb0), self.from_rgb1(srgb1)]

        # linspace
        ls = np.linspace(cs[0], cs[1], endpoint=True, num=n, axis=0)

        # back to srgb
        srgb = self.to_rgb1(ls.T).T

        srgb[srgb < 0] = 0.0
        srgb[srgb > 1] = 1.0
        return srgb

    def show_srgb255_gradient(self, srgb0, srgb1, n=256):
        srgb0 = np.asarray(srgb0)
        srgb1 = np.asarray(srgb1)
        self.show_srgb1_gradient(srgb0 / 255, srgb1 / 255, n)

    def plot_srgb255_gradient(self, srgb0, srgb1, n=256):
        srgb0 = np.asarray(srgb0)
        srgb1 = np.asarray(srgb1)
        self.plot_srgb1_gradient(srgb0 / 255, srgb1 / 255, n)

    def get_srgb255_gradient(self, srgb0, srgb1, n):
        srgb0 = np.asarray(srgb0)
        srgb1 = np.asarray(srgb1)
        return self.get_srgb1_gradient(srgb0 / 255, srgb1 / 255, n) * 255

    def show_primary_srgb_gradients(self, *args, **kwargs):
        self.plot_primary_srgb_gradients(*args, **kwargs)
        plt.show()
        plt.close()

    def plot_primary_srgb_gradients(self, n=256):
        pairs = [
            [([1, 1, 1], [1, 0, 0]), ([1, 0, 0], [0, 1, 0])],
            [([1, 1, 1], [0, 1, 0]), ([0, 1, 0], [0, 0, 1])],
            [([1, 1, 1], [0, 0, 1]), ([0, 0, 1], [1, 0, 0])],
            [([0, 0, 0], [1, 0, 0]), ([1, 0, 0], [0, 1, 1])],
            [([0, 0, 0], [0, 1, 0]), ([0, 1, 0], [1, 0, 1])],
            [([0, 0, 0], [0, 0, 1]), ([0, 0, 1], [1, 1, 0])],
        ]
        fig, axes = plt.subplots(len(pairs), 2)
        for i in range(len(pairs)):
            for j in range(2):
                pair = pairs[i][j]
                ax = axes[i][j]
                srgb = self.get_srgb1_gradient(pair[0], pair[1], n=n)

                cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", srgb, n)

                gradient = np.linspace(0.0, 1.0, n)
                gradient = np.vstack((gradient, gradient))
                ax.imshow(gradient, aspect="auto", cmap=cmap)
                ax.axis("off")
        fig.suptitle(f"primary SRGB gradients in {self.name}")


def _xyy_to_xyz100(xyy):
    x, y, Y = xyy
    return np.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100
