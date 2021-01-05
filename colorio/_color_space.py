import os

import matplotlib
import matplotlib.pyplot as plt
import numpy
import yaml

from ._hdr import HdrLinear
from ._srgb import SrgbLinear
from ._tools import get_mono_outline_xy, get_munsell_data, spectrum_to_xyz100
from .illuminants import whitepoints_cie1931
from .observers import cie_1931_2


class ColorSpace:
    def save_visible_gamut(self, observer, illuminant, filename, cut_000=False):
        import meshio
        from scipy.spatial import ConvexHull

        lmbda, illu = illuminant
        values = []

        # Iterate over every possible illuminant input and store it in values
        n = len(lmbda)
        values = numpy.empty((n * (n - 1) + 2, 3))
        k = 0

        # No light
        data = numpy.zeros(n)
        values[k] = spectrum_to_xyz100((lmbda, illu * data), observer=observer)
        k += 1
        # frequency blocks
        for width in range(1, n):
            data = numpy.zeros(n)
            data[:width] = 1.0
            for _ in range(n):
                values[k] = spectrum_to_xyz100((lmbda, illu * data), observer=observer)
                k += 1
                data = numpy.roll(data, shift=1)
        # Full illuminant
        data = numpy.ones(len(lmbda))
        values[k] = spectrum_to_xyz100((lmbda, illu * data), observer=observer)
        k += 1

        # scale the values such that the Y-coordinate of the white point (last entry)
        # has value 100.
        values *= 100 / values[-1][1]

        cells = ConvexHull(values).simplices

        if cut_000:
            values = values[1:]
            cells = cells[~numpy.any(cells == 0, axis=1)]
            cells -= 1

        pts = self.from_xyz100(values.T).T

        meshio.write_points_cells(filename, pts, cells={"triangle": cells})

    def save_srgb_gamut(self, filename, n=50, cut_000=False):
        import meshio
        import meshzoo

        points, cells = meshzoo.cube(nx=n, ny=n, nz=n)

        if cut_000:
            # cut off [0, 0, 0] to avoid division by 0 in the xyz conversion
            points = points[1:]
            cells = cells[~numpy.any(cells == 0, axis=1)]
            cells -= 1

        srgb_linear = SrgbLinear()
        pts = self.from_xyz100(srgb_linear.to_xyz100(points.T)).T
        assert pts.shape[1] == 3
        rgb = srgb_linear.to_srgb1(points)
        meshio.write_points_cells(
            filename, pts, {"tetra": cells}, point_data={"srgb": rgb}
        )

    def save_hdr_gamut(self, filename, n=50, cut_000=False):
        import meshio
        import meshzoo

        points, cells = meshzoo.cube(nx=n, ny=n, nz=n)

        if cut_000:
            # cut off [0, 0, 0] to avoid division by 0 in the xyz conversion
            points = points[1:]
            cells = cells[~numpy.any(cells == 0, axis=1)]
            cells -= 1

        hdr = HdrLinear()
        pts = self.from_xyz100(hdr.to_xyz100(points.T)).T
        rgb = hdr.to_hdr1(points)
        meshio.write_points_cells(
            filename, pts, {"tetra": cells}, point_data={"hdr-rgb": rgb}
        )

    def save_cone_gamut(self, filename, observer, max_Y):
        import meshio
        import pygmsh

        with pygmsh.geo.Geometry() as geom:
            max_stepsize = 4.0e-2
            xy, _ = get_mono_outline_xy(observer, max_stepsize=max_stepsize)

            # append third component
            xy = numpy.column_stack([xy, numpy.full(xy.shape[0], 1.0e-5)])

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

    def show_luo_rigg(self, *args, **kwargs):
        plt.figure()
        self.plot_luo_rigg(*args, **kwargs)
        plt.show()
        plt.close()

    def save_luo_rigg(self, filename, *args, **kwargs):
        plt.figure()
        self.plot_luo_rigg(*args, **kwargs)
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        plt.close()

    def plot_luo_rigg(
        self, level, outline_prec=1.0e-2, plot_srgb_gamut=True, ellipse_scaling=2.0
    ):
        # M. R. Luo, B. Rigg,
        # Chromaticity Discrimination Ellipses for Surface Colours,
        # Color Research and Application, Volume 11, Issue 1, Spring 1986, Pages 25-42,
        # <https://doi.org/10.1002/col.5080110107>.
        dir_path = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dir_path, "data/luo-rigg/luo-rigg.yaml")) as f:
            data = yaml.safe_load(f)
        #
        xy_centers = []
        xy_offsets = []
        # collect the ellipse centers and offsets
        # Use four offset points of each ellipse, one could take more
        alpha = 2 * numpy.pi * numpy.linspace(0.0, 1.0, 16, endpoint=False)
        pts = numpy.array([numpy.cos(alpha), numpy.sin(alpha)])
        for data_set in data.values():
            # The set factor is the mean of the R values
            # set_factor = sum([dat[-1] for dat in data_set.values()]) / len(data_set)
            for dat in data_set.values():
                x, y, Y, a, a_div_b, theta_deg, _ = dat
                theta = theta_deg * 2 * numpy.pi / 360
                a /= 1.0e4
                a *= (Y / 30) ** 0.2
                b = a / a_div_b
                # plot the ellipse
                xy_centers.append([x, y])
                J = numpy.array(
                    [
                        [+a * numpy.cos(theta), -b * numpy.sin(theta)],
                        [+a * numpy.sin(theta), +b * numpy.cos(theta)],
                    ]
                )
                xy_offsets.append(numpy.dot(J, pts))

        self._plot_ellipses(
            xy_centers,
            xy_offsets,
            level,
            outline_prec=outline_prec,
            plot_srgb_gamut=plot_srgb_gamut,
            ellipse_scaling=ellipse_scaling,
        )

    def show_munsell(self, V):
        plt.figure()
        self.plot_munsell(V)
        plt.show()
        plt.close()

    def save_munsell(self, filename, V):
        plt.figure()
        self.plot_munsell(V)
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        plt.close()

    def plot_munsell(self, V):
        _, v, _, xyy = get_munsell_data()

        # pick the data from the given munsell level
        xyy = xyy[:, v == V]

        x, y, Y = xyy
        xyz100 = numpy.array([Y / y * x, Y, Y / y * (1 - x - y)])
        vals = self.from_xyz100(xyz100)

        srgb = SrgbLinear()
        rgb = srgb.from_xyz100(xyz100)
        is_legal_srgb = numpy.all((0 <= rgb) & (rgb <= 1), axis=0)

        idx = [0, 1, 2]
        k1, k2 = idx[: self.k0] + idx[self.k0 + 1 :]

        # plot the ones that cannot be represented in SRGB
        plt.plot(
            vals[k1, ~is_legal_srgb],
            vals[k2, ~is_legal_srgb],
            "o",
            color="white",
            markeredgecolor="black",
        )
        # plot the srgb dots
        for val, rgb_ in zip(vals[:, is_legal_srgb].T, rgb[:, is_legal_srgb].T):
            plt.plot(val[k1], val[k2], "o", color=srgb.to_srgb1(rgb_))

        plt.grid()
        plt.title(f"V={V}")
        plt.xlabel(self.labels[k1])
        plt.ylabel(self.labels[k2])
        plt.axis("equal")
