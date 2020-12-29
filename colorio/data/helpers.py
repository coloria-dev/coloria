import numpy
import matplotlib.pyplot as plt

from .._srgb import SrgbLinear


def _plot_color_constancy_data(
    data_xyz100, wp_xyz100, colorspace, approximate_colors_in_srgb=False
):
    # k0 is the coordinate that corresponds to "lightness"
    k0 = colorspace.k0

    k1, k2 = [k for k in [0, 1, 2] if k != k0]

    wp = colorspace.from_xyz100(wp_xyz100)[[k1, k2]]
    srgb = SrgbLinear()
    for xyz in data_xyz100:
        d = colorspace.from_xyz100(xyz)[[k1, k2]]

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

    plt.xlabel(colorspace.labels[k1])
    plt.ylabel(colorspace.labels[k2])
    plt.grid()
    plt.axis("equal")
