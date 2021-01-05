import matplotlib
import matplotlib.pyplot as plt
import numpy

from ._srgb import SrgbLinear


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


def _plot_srgb_slice(cs, lightness):
    import meshzoo

    # Get all RGB values that sum up to 1.
    bary, triangles = meshzoo.triangle(n=50)
    corners = numpy.array([[0, 0], [1, 0], [0, 1]]).T
    srgb_vals = numpy.dot(corners, bary).T
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
        self_val_min = cs.from_xyz100(xyz100)[cs.k0]
        if self_val_min > lightness:
            mask[k] = False
            continue

        alpha_max = 1.0 / numpy.max(val)

        xyz100 = srgb_linear.to_xyz100(val * alpha_max)
        self_val_max = cs.from_xyz100(xyz100)[cs.k0]
        if self_val_max < lightness:
            mask[k] = False
            continue

        while True:
            alpha = (alpha_max + alpha_min) / 2
            srgb_linear_vals[k] = val * alpha
            xyz100 = srgb_linear.to_xyz100(srgb_linear_vals[k])
            self_val = cs.from_xyz100(xyz100)
            if abs(self_val[cs.k0] - lightness) < tol:
                break
            elif self_val[cs.k0] < lightness:
                alpha_min = alpha
            else:
                assert self_val[cs.k0] > lightness
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
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("gamut", rgb, N=len(rgb))

    k1, k2 = [k for k in [0, 1, 2] if k != cs.k0]

    plt.tripcolor(
        self_vals[:, k1],
        self_vals[:, k2],
        triangles,
        z,
        shading="gouraud",
        cmap=cmap,
    )
    # plt.triplot(self_vals[:, k1], self_vals[:, k2], triangles=triangles)


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
        plot_srgb_gamut=plot_srgb_gamut,
        fill_color=visible_gamut_fill_color,
    )

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
