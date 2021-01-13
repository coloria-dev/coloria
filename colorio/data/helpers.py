import matplotlib.pyplot as plt
import numpy

from ..cs import SrgbLinear


def _compute_straight_line_residuals(cs, wp, d):
    """Compute the TLS residuals for each of the arms."""
    # remove the row corresponding to lightness
    idx = [0, 1, 2]
    idx.pop(cs.k0)
    wp_cs = cs.from_xyz100(wp)[idx]
    s2 = []
    for dd in d:
        vals = cs.from_xyz100(dd)[idx]
        # move values such that whitepoint is in the origin
        vals = (vals.T - wp_cs).T
        # could also be computed explicitly
        s_max, s_min = numpy.linalg.svd(vals, compute_uv=False)
        s2.append(s_min / s_max)
        # plt.plot(vals[0], vals[1], "x")
        # plt.gca().set_aspect("equal")
        # plt.show()
    return numpy.array(s2)


def _plot_hue_linearity_data(
    data_xyz100, wp_xyz100, colorspace, approximate_colors_in_srgb=False
):
    # k0 is the coordinate that corresponds to "lightness"
    k0 = colorspace.k0
    k1, k2 = [k for k in [0, 1, 2] if k != k0]

    wp = colorspace.from_xyz100(wp_xyz100)[[k1, k2]]
    srgb = SrgbLinear()
    for xyz in data_xyz100:
        d = colorspace.from_xyz100(xyz)[[k1, k2]]

        # get the eigenvector corresponding to the larger eigenvalue
        d_wp = (d.T - wp).T
        vals, vecs = numpy.linalg.eigh(d_wp @ d_wp.T)
        v = vecs[:, 0] if vals[0] > vals[1] else vecs[:, 1]

        if numpy.dot(v, numpy.average(d, axis=1)) < 0:
            v = -v

        length = numpy.sqrt(numpy.max(numpy.einsum("ij,ij->i", d.T - wp, d.T - wp)))
        end_point = wp + length * v
        plt.plot([wp[0], end_point[0]], [wp[1], end_point[1]], "-", color="0.5")

        for dd, rgb in zip(d.T, srgb.from_xyz100(xyz).T):
            if approximate_colors_in_srgb:
                is_legal_srgb = True
                rgb[rgb > 1] = 1
                rgb[rgb < 0] = 0
            else:
                is_legal_srgb = numpy.all(rgb >= 0) and numpy.all(rgb <= 1)
            col = srgb.to_rgb1(rgb) if is_legal_srgb else "white"
            ecol = srgb.to_rgb1(rgb) if is_legal_srgb else "black"
            plt.plot(dd[0], dd[1], "o", color=col, markeredgecolor=ecol)

    plt.xlabel(colorspace.labels[k1])
    plt.ylabel(colorspace.labels[k2])
    plt.axis("equal")

    # plt.grid()
    plt.grid(False)
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)


def _compute_ellipse_residual(cs, xyy100_centers, xyy100_points):
    distances = []
    for xyy100_center, xyy100_pts in zip(xyy100_centers, xyy100_points):
        # append Y
        xyz100_center = _xyy100_to_xyz100(xyy100_center)
        assert numpy.all(xyz100_center >= 0.0)
        xyz100_pts = _xyy100_to_xyz100(xyy100_pts)
        assert numpy.all(xyz100_pts >= 0.0)
        # plt.plot(xyz_center[0], xyz_center[2], "x")
        # plt.plot(xyz_ellips[0], xyz_ellips[2], "o")
        # plt.show()
        cs_center = cs.from_xyz100(xyz100_center)
        cs_ellips = cs.from_xyz100(xyz100_pts)
        # plt.plot(cs_center[0], cs_center[1], "x")
        # plt.plot(cs_ellips[0], cs_ellips[1], "o")
        # plt.show()
        #
        # compute distances to ellipse center
        diff = (cs_center - cs_ellips.T).T
        distances.append(numpy.sqrt(diff[0] ** 2 + diff[1] ** 2))

    distances = numpy.concatenate(distances)
    alpha = numpy.average(distances)
    return numpy.sqrt(numpy.sum((alpha - distances) ** 2) / numpy.sum(distances ** 2))


def _xyy100_to_xyz100(xyy):
    x, y, Y = xyy
    return numpy.array([Y / y * x, Y, Y / y * (1 - x - y)])


def _plot_ellipses(xyy100_centers, xyy100_points, cs, ellipse_scaling=1.0):
    from matplotlib.patches import Ellipse
    from scipy.optimize import leastsq

    for center, points in zip(xyy100_centers, xyy100_points):
        # cut off the irrelevant index
        cs_center = cs.from_xyz100(_xyy100_to_xyz100(center))
        cs_points = cs.from_xyz100(_xyy100_to_xyz100(points))

        # project out lightness component
        tcenter = numpy.delete(cs_center, cs.k0)
        tvals = numpy.delete(cs_points, cs.k0, axis=0)

        # Given these new transformed vals, find the ellipse that best fits those
        # points
        X = (tvals.T - tcenter).T

        def f_ellipse(a_b_theta):
            a, b, theta = a_b_theta
            sin_t = numpy.sin(theta)
            cos_t = numpy.cos(theta)
            return (
                +(a ** 2) * (X[0] * cos_t + X[1] * sin_t) ** 2
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
            angle=theta / numpy.pi * 180,
            # label=label,
        )
        plt.gca().add_artist(e)
        e.set_alpha(0.5)
        e.set_facecolor("k")

        # plt.plot(tvals[:, 0], tvals[:, 1], "xk")
    plt.gca().set_aspect("equal")
    labels = cs.labels[: cs.k0] + cs.labels[cs.k0 + 1 :]
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])

    # mpl doesn't update axis limits when adding artists,
    # <https://github.com/matplotlib/matplotlib/issues/19290>.
    # Handle it manually.
    tcenters = []
    for center, points in zip(xyy100_centers, xyy100_points):
        cs_center = cs.from_xyz100(_xyy100_to_xyz100(center))
        tcenters.append(numpy.delete(cs_center, cs.k0))
    tcenters = numpy.asarray(tcenters).T
    xmin = numpy.min(tcenters[0])
    xmax = numpy.max(tcenters[0])
    ymin = numpy.min(tcenters[1])
    ymax = numpy.max(tcenters[1])
    width = xmax - xmin
    height = ymax - ymin
    plt.xlim(xmin - 0.2 * width, xmax + 0.2 * width)
    plt.ylim(ymin - 0.2 * height, ymax + 0.2 * height)
