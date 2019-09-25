import os

import matplotlib
import matplotlib.pyplot as plt
import numpy
import yaml
from matplotlib.patches import Ellipse
from scipy.optimize import leastsq
from scipy.spatial import ConvexHull

from ._hdr import HdrLinear
from ._srgb import SrgbLinear
from ._tools import get_mono_outline_xy, get_munsell_data, spectrum_to_xyz100
from .illuminants import whitepoints_cie1931
from .observers import cie_1931_2


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

    def save_cone_gamut(self, filename, observer, max_Y):
        import meshio
        import pygmsh

        geom = pygmsh.built_in.Geometry()

        max_stepsize = 4.0e-2
        xy, _ = get_mono_outline_xy(observer, max_stepsize=max_stepsize)

        # append third component
        xy = numpy.column_stack([xy, numpy.full(xy.shape[0], 1.0e-5)])

        # Draw a cross.
        poly = geom.add_polygon(xy, lcar=max_stepsize)

        axis = [0, 0, max_Y]

        geom.extrude(poly, translation_axis=axis, point_on_axis=[0, 0, 0])

        mesh = pygmsh.generate_mesh(geom, verbose=False)
        # meshio.write(filename, mesh)

        pts = self.from_xyz100(_xyy_to_xyz100(mesh.points.T)).T
        meshio.write_points_cells(filename, pts, {"tetra": mesh.cells["tetra"]})
        return

    def show_macadam(self, *args, **kwargs):
        self.plot_macadam(*args, **kwargs)
        plt.show()
        return

    def save_macadam(self, filename, *args, **kwargs):
        self.plot_macadam(*args, **kwargs)
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        return

    def plot_macadam(
        self, k0, level, outline_prec=1.0e-2, plot_srgb_gamut=True, ellipse_scaling=10.0
    ):
        # Extract ellipse centers and offsets from MacAdams data
        dir_path = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dir_path, "data/macadam1942/table3.yaml")) as f:
            data = yaml.safe_load(f)
        #
        # collect the ellipse centers and offsets
        xy_centers = []
        xy_offsets = []
        for datak in data:
            # collect ellipse points
            _, _, _, _, delta_y_delta_x, delta_s = numpy.array(datak["data"]).T
            offset = (
                numpy.array([numpy.ones(delta_y_delta_x.shape[0]), delta_y_delta_x])
                / numpy.sqrt(1 + delta_y_delta_x ** 2)
                * delta_s
            )
            if offset.shape[1] < 2:
                continue
            xy_centers.append([datak["x"], datak["y"]])
            xy_offsets.append(numpy.column_stack([+offset, -offset]))

        self._plot_ellipses(
            xy_centers,
            xy_offsets,
            k0,
            level,
            outline_prec=outline_prec,
            plot_srgb_gamut=plot_srgb_gamut,
            ellipse_scaling=ellipse_scaling,
        )
        return

    def show_luo_rigg(self, *args, **kwargs):
        self.plot_luo_rigg(*args, **kwargs)
        plt.show()
        return

    def save_luo_rigg(self, filename, *args, **kwargs):
        self.plot_luo_rigg(*args, **kwargs)
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        return

    def plot_luo_rigg(
        self, k0, level, outline_prec=1.0e-2, plot_srgb_gamut=True, ellipse_scaling=2.0
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
            k0,
            level,
            outline_prec=outline_prec,
            plot_srgb_gamut=plot_srgb_gamut,
            ellipse_scaling=ellipse_scaling,
        )
        return

    def _plot_ellipses(
        self,
        xy_centers,
        xy_offsets,
        k0,
        level,
        outline_prec=1.0e-2,
        plot_srgb_gamut=True,
        ellipse_scaling=10.0,
    ):
        self.plot_visible_slice(
            k0, level, outline_prec=outline_prec, plot_srgb_gamut=plot_srgb_gamut
        )

        for center, offset in zip(xy_centers, xy_offsets):
            # get all the approximate ellipse points in xy space
            xy_ellipse = (center + offset.T).T

            tcenter = self._bisect(center, k0, level)
            tvals = numpy.array([self._bisect(xy, k0, level) for xy in xy_ellipse.T])

            # cut off the irrelevant index
            idx = [0, 1, 2]
            k1, k2 = idx[:k0] + idx[k0 + 1 :]
            tcenter = tcenter[[k1, k2]]
            tvals = tvals[:, [k1, k2]]

            # Given these new transformed vals, find the ellipse that best fits those
            # points
            X = (tvals - tcenter).T

            def f_ellipse(a_b_theta):
                a, b, theta = a_b_theta
                sin_t = numpy.sin(theta)
                cos_t = numpy.cos(theta)
                return (
                    +a ** 2 * (X[0] * cos_t + X[1] * sin_t) ** 2
                    + b ** 2 * (X[0] * sin_t - X[1] * cos_t) ** 2
                    - 1.0
                )

            def jac(a_b_theta):
                a, b, theta = a_b_theta
                x0sin = X[0] * numpy.sin(theta)
                x0cos = X[0] * numpy.cos(theta)
                x1sin = X[1] * numpy.sin(theta)
                x1cos = X[1] * numpy.cos(theta)
                return numpy.array(
                    [
                        +2 * a * (x0cos + x1sin) ** 2,
                        +2 * b * (x0sin - x1cos) ** 2,
                        +a ** 2 * 2 * (x0cos + x1sin) * (-x0sin + x1cos)
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
                angle=theta / numpy.pi * 180,
                # label=label,
            )
            plt.gca().add_artist(e)
            e.set_alpha(0.5)
            e.set_facecolor("k")

            # plt.plot(tvals[:, 0], tvals[:, 1], "xk")
        return

    def show_visible_slice(self, *args, **kwargs):
        self.plot_visible_slice(*args, **kwargs)
        plt.show()
        return

    def save_visible_slice(self, filename, *args, **kwargs):
        self.plot_visible_slice(*args, **kwargs)
        plt.savefig(filename, transparent=True, bbox_inches="tight")
        return

    def plot_visible_slice(self, k0, level, outline_prec=1.0e-2, plot_srgb_gamut=True):
        # first plot the monochromatic outline
        mono_xy, conn_xy = get_mono_outline_xy(
            observer=cie_1931_2(), max_stepsize=outline_prec
        )

        mono_vals = numpy.array([self._bisect(xy, k0, level) for xy in mono_xy])
        conn_vals = numpy.array([self._bisect(xy, k0, level) for xy in conn_xy])

        idx = [0, 1, 2]
        k1, k2 = idx[:k0] + idx[k0 + 1 :]
        plt.plot(mono_vals[:, k1], mono_vals[:, k2], "-", color="k")
        plt.plot(conn_vals[:, k1], conn_vals[:, k2], ":", color="k")
        #
        xyz = numpy.vstack([mono_vals, conn_vals[1:]])
        plt.fill(xyz[:, k1], xyz[:, k2], facecolor="0.8", zorder=0)

        if plot_srgb_gamut:
            self._plot_srgb_gamut(k0, level)

        plt.axis("equal")
        plt.xlabel(self.labels[k1])
        plt.ylabel(self.labels[k2])
        plt.title(f"{self.labels[k0]} = {level}")
        return

    def _plot_srgb_gamut(self, k0, level, bright=False):
        import meshzoo

        # Get all RGB values that sum up to 1.
        srgb_vals, triangles = meshzoo.triangle(n=50, corners=[[0, 0], [1, 0], [0, 1]])
        srgb_vals = numpy.column_stack([srgb_vals, 1.0 - numpy.sum(srgb_vals, axis=1)])

        # matplotlib is sensitive when it comes to srgb values, so take good care here
        assert numpy.all(srgb_vals > -1.0e-10)
        srgb_vals[srgb_vals < 0.0] = 0.0

        # Use bisection to
        srgb_linear = SrgbLinear()
        tol = 1.0e-5
        # Use zeros() instead of empty() here to avoid invalid values when setting up
        # the cmap below.
        self_vals = numpy.zeros((srgb_vals.shape[0], 3))
        srgb_linear_vals = numpy.zeros((srgb_vals.shape[0], 3))
        mask = numpy.ones(srgb_vals.shape[0], dtype=bool)
        for k, val in enumerate(srgb_vals):
            alpha_min = 0.0
            xyz100 = srgb_linear.to_xyz100(val * alpha_min)
            self_val_min = self.from_xyz100(xyz100)
            if self_val_min[k0] > level:
                mask[k] = False
                continue

            alpha_max = 1.0 / numpy.max(val)

            xyz100 = srgb_linear.to_xyz100(val * alpha_max)
            self_val_max = self.from_xyz100(xyz100)
            if self_val_max[k0] < level:
                mask[k] = False
                continue

            while True:
                alpha = (alpha_max + alpha_min) / 2
                srgb_linear_vals[k] = val * alpha
                xyz100 = srgb_linear.to_xyz100(srgb_linear_vals[k])
                self_val = self.from_xyz100(xyz100)
                if abs(self_val[k0] - level) < tol:
                    break
                elif self_val[k0] < level:
                    alpha_min = alpha
                else:
                    assert self_val[k0] > level
                    alpha_max = alpha
            self_vals[k] = self_val

        # Remove all triangles which have masked corner points
        tri_mask = numpy.all(mask[triangles], axis=1)
        if ~numpy.any(tri_mask):
            return
        triangles = triangles[tri_mask]

        # Unfortunately, one cannot use tripcolors with explicit RGB specification (see
        # <https://github.com/matplotlib/matplotlib/issues/10265>). As a workaround,
        # associate range(n) data with the points and create a colormap that associates
        # the integer values with the respective RGBs.
        z = numpy.arange(srgb_vals.shape[0])
        rgb = srgb_linear.to_srgb1(srgb_linear_vals)
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "gamut", rgb, N=len(rgb)
        )

        idx = [0, 1, 2]
        k1, k2 = idx[:k0] + idx[k0 + 1 :]

        plt.tripcolor(
            self_vals[:, k1],
            self_vals[:, k2],
            triangles,
            z,
            shading="gouraud",
            cmap=cmap,
        )
        # plt.triplot(self_vals[:, k1], self_vals[:, k2], triangles=triangles)
        return

    def _bisect(self, xy, k0, level):
        tol = 1.0e-5
        x, y = xy
        # Use bisection to find a matching Y value that projects the xy into the
        # given level.
        min_Y = 0.0
        xyz100 = numpy.array([min_Y / y * x, min_Y, min_Y / y * (1 - x - y)]) * 100
        min_val = self.from_xyz100(xyz100)[k0]
        assert min_val <= level

        # search for an appropriate max_Y to start with
        max_Y = 1.0
        while True:
            xyz100 = numpy.array([max_Y / y * x, max_Y, max_Y / y * (1 - x - y)]) * 100
            max_val = self.from_xyz100(xyz100)[k0]
            if max_val >= level:
                break
            max_Y *= 2

        while True:
            Y = (max_Y + min_Y) / 2
            xyz100 = numpy.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100
            val = self.from_xyz100(xyz100)
            if abs(val[k0] - level) < tol:
                break
            elif val[k0] > level:
                max_Y = Y
            else:
                assert val[k0] < level
                min_Y = Y

        return val

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

        wp = numpy.array(data["white point"])

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

        wp = numpy.array(whitepoints_cie1931["C"])
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


def _plot_color_constancy_data(
    data_xyz100, wp_xyz100, colorspace, approximate_colors_in_srgb=False
):
    wp = colorspace.from_xyz100(wp_xyz100)[1:]
    srgb = SrgbLinear()
    for xyz in data_xyz100:
        d = colorspace.from_xyz100(xyz)[1:]

        # There are numerous possibilities of defining the "best" approximating line for
        # a bunch of points (x_i, y_i). For example, one could try and minimize the
        # expression
        #    sum_i (-numpy.sin(theta) * x_i + numpy.cos(theta) * y_i) ** 2
        # over theta, which means to minimize the orthogonal component of (x_i, y_i) to
        # (cos(theta), sin(theta)).
        #
        # A more simple and effective approach is to use the average of all points,
        #    theta = arctan(sum(y_i) / sum(x_i)).
        # This also fits in nicely with minimization problems which move around the
        # points to minimize the difference from the average,
        #
        #    sum_j (y_j / x_j - bar{y} / bar{x}) ** 2 -> min,
        #    sum_j (y_j bar{x} - x_j bar{y}) ** 2 -> min.
        #
        # Plot it from wp to the outmost point
        avg = numpy.sum(d, axis=1) / d.shape[1]
        length = numpy.sqrt(numpy.max(numpy.einsum("ij,ij->i", d.T - wp, d.T - wp)))
        end_point = wp + length * (avg - wp) / numpy.sqrt(numpy.sum((avg - wp) ** 2))
        plt.plot([wp[0], end_point[0]], [wp[1], end_point[1]], "-", color="0.5")

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

    plt.xlabel(colorspace.labels[1])
    plt.ylabel(colorspace.labels[2])
    plt.grid()
    plt.axis("equal")
    return


def _xyy_from_xyz100(xyz):
    sum_xyz = numpy.sum(xyz, axis=0)
    x = xyz[0]
    y = xyz[1]
    return numpy.array([x / sum_xyz, y / sum_xyz, y / 100])


def _xyy_to_xyz100(xyy):
    x, y, Y = xyy
    return numpy.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100
