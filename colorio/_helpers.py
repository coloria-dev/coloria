import matplotlib.pyplot as plt
import numpy


def _find_Y(cs, xy, level, tol=1.0e-5):
    """Use bisection to find a matching Y value that projects the xy into the given
    level.
    """
    x, y = xy
    min_Y = 0.0
    xyz100 = numpy.array([min_Y / y * x, min_Y, min_Y / y * (1 - x - y)]) * 100
    min_val = cs.from_xyz100(xyz100)[cs.k0]
    assert min_val <= level

    # search for an appropriate max_Y to start with
    max_Y = 1.0
    while True:
        xyz100 = numpy.array([max_Y / y * x, max_Y, max_Y / y * (1 - x - y)]) * 100
        max_val = cs.from_xyz100(xyz100)[cs.k0]
        if max_val >= level:
            break
        max_Y *= 2

    while True:
        Y = (max_Y + min_Y) / 2
        xyz100 = numpy.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100
        val = cs.from_xyz100(xyz100)
        if abs(val[cs.k0] - level) < tol:
            break
        elif val[cs.k0] > level:
            max_Y = Y
        else:
            assert val[cs.k0] < level
            min_Y = Y

    return val


def _plot_ellipses(
    xy_centers,
    xy_offsets,
    cs,
    lightness,
    outline_prec=1.0e-2,
    plot_srgb_gamut=True,
    ellipse_scaling=10.0,
    visible_gamut_fill_color="0.8",
):
    from matplotlib.patches import Ellipse
    from scipy.optimize import leastsq

    cs.plot_visible_slice(
        lightness,
        outline_prec=outline_prec,
        fill_color=visible_gamut_fill_color,
    )
    if plot_srgb_gamut:
        cs.plot_srgb_slice(lightness)

    for center, offset in zip(xy_centers, xy_offsets):
        # get all the approximate ellipse points in xy space
        xy_ellipse = (center + offset.T).T

        tcenter = _find_Y(cs, center, lightness)
        tvals = numpy.array([_find_Y(cs, xy, lightness) for xy in xy_ellipse.T])

        # cut off the irrelevant index
        k1, k2 = [k for k in [0, 1, 2] if k != cs.k0]
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
