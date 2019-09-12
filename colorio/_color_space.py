import os

import matplotlib.pyplot as plt
import numpy
import yaml
from scipy.optimize import leastsq
from scipy.spatial import ConvexHull

from ._hdr import HdrLinear
from ._srgb import SrgbLinear
from ._tools import get_munsell_data, spectrum_to_xyz100
from .illuminants import whitepoints_cie1931


class ColorSpace:
    def __init__(self):
        return

    def save_visible_gamut(self, observer, illuminant, filename, cut_000=False):
        import meshio

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
        return

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
        return

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
        return

    def show_ebner_fairchild(self):
        self.plot_ebner_fairchild()
        plt.show()
        return

    def save_ebner_fairchild(self, filename):
        self.plot_ebner_fairchild()
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        return

    def plot_ebner_fairchild(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dir_path, "data/ebner_fairchild.yaml")) as f:
            data = yaml.safe_load(f)

        wp = self.from_xyz100(numpy.array(data["white point"]))[1:]

        d = [
            numpy.column_stack([dat["reference xyz"], numpy.array(dat["same"]).T])
            for dat in data["data"]
        ]

        _plot_color_constancy_data(d, wp, self)
        return

    def show_hung_berns(self):
        self.plot_hung_berns()
        plt.show()
        return

    def save_hung_berns(self, filename):
        self.plot_hung_berns()
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        return

    def plot_hung_berns(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dir_path, "data/hung-berns/table3.yaml")) as f:
            data = yaml.safe_load(f)

        wp = self.from_xyz100(numpy.array(whitepoints_cie1931["C"]))[1:]

        d = [numpy.array(list(color.values())).T for color in data.values()]

        _plot_color_constancy_data(d, wp, self)
        return

    def show_xiao(self):
        self.plot_xiao()
        plt.show()
        return

    def save_xiao(self, filename):
        self.plot_xiao()
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        return

    def plot_xiao(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))

        filenames = [
            "unique_blue.yaml",
            "unique_green.yaml",
            "unique_red.yaml",
            "unique_yellow.yaml",
        ]

        data = []
        for filename in filenames:
            with open(os.path.join(dir_path, "data", "xiao", filename)) as f:
                dat = numpy.array(yaml.safe_load(f))
            # average over all observers and sessions
            data.append(numpy.sum(dat, axis=(0, 1)) / numpy.prod(dat.shape[:2]))

        data = numpy.array(data)

        # Use Xiao's 'neutral gray' as white point.
        with open(os.path.join(dir_path, "data/xiao/neutral_gray.yaml")) as f:
            ng_data = numpy.array(yaml.safe_load(f))

        ng = numpy.sum(ng_data, axis=0) / numpy.prod(ng_data.shape[:1])
        ng_cs = self.from_xyz100(ng)[1:]

        data = numpy.moveaxis(data, 1, 2)

        _plot_color_constancy_data(data, ng_cs, self)
        return

    def show_munsell(self, V):
        self.plot_munsell(V)
        plt.show()
        return

    def save_munsell(self, filename, V):
        self.plot_munsell(V)
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        return

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

        # plot the ones that cannot be represented in SRGB
        plt.plot(
            vals[1, ~is_legal_srgb],
            vals[2, ~is_legal_srgb],
            "o",
            color="white",
            markeredgecolor="black",
        )
        # plot the srgb dots
        for val, rgb_ in zip(vals[:, is_legal_srgb].T, rgb[:, is_legal_srgb].T):
            plt.plot(val[1], val[2], "o", color=srgb.to_srgb1(rgb_))

        plt.axis("equal")
        return


def _plot_color_constancy_data(data, wp, colorspace, approximate_colors_in_srgb=False):
    srgb = SrgbLinear()
    for xyz in data:
        d = colorspace.from_xyz100(xyz)[1:]

        # Find best fit line through all points
        def f(theta, D=d):
            return +numpy.sin(theta) * (D[0] - wp[0]) + numpy.cos(theta) * (
                D[1] - wp[1]
            )

        def jac(theta, D=d):
            return +numpy.cos(theta) * (D[0] - wp[0]) - numpy.sin(theta) * (
                D[1] - wp[1]
            )

        # out = least_squares(f, 0.0)
        # out = leastsq(f, 0.0, full_output=True)

        out, _ = leastsq(f, 0.0, Dfun=jac)
        # We have to take the first element here, see
        # <https://github.com/scipy/scipy/issues/8532>.
        theta = out[0]

        # Plot it from wp to the outmost point
        length = numpy.sqrt(numpy.max(numpy.einsum("ij,ij->i", (d.T - wp), (d.T - wp))))
        # The solution theta can be rotated by pi and still give the same
        # result. Find out on which side all the points are sitting and plot
        # the line accordingly.
        ex = length * numpy.array([numpy.cos(theta), -numpy.sin(theta)])

        end_point = wp + ex
        ep_d = numpy.linalg.norm(end_point - d[:, -1])
        ep_wp = numpy.linalg.norm(end_point - wp)
        if ep_d > ep_wp:
            end_point = wp - ex
        plt.plot([wp[0], end_point[0]], [wp[1], end_point[1]], "-", color="0.5")

        # Deliberatly only handle the two last components, e.g., a* b* from
        # L*a*b*. They typically indicate the chroma.
        for dd, rgb in zip(d.T, srgb.from_xyz100(xyz).T):
            if approximate_colors_in_srgb:
                is_legal_srgb = True
                rgb[rgb > 1] = 1
                rgb[rgb < 0] = 0
            else:
                is_legal_srgb = numpy.all(rgb >= 0) and numpy.all(rgb <= 1)
            col = srgb.to_srgb1(rgb) if is_legal_srgb else "white"
            ecol = srgb.to_srgb1(rgb) if is_legal_srgb else "black"
            plt.plot(dd[0], dd[1], "o", color=col, markeredgecolor=ecol)

    plt.axis("equal")
    return
