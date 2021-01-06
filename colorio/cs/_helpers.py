import matplotlib
import matplotlib.pyplot as plt
import meshzoo
import numpy


def _plot_rgb_triangle(xy_to_2d, bright=True):
    # plot sRGB triangle
    # discretization points
    n = 50

    # Get all RGB values that sum up to 1.
    rgb_linear, _ = meshzoo.triangle(n)
    if bright:
        # For the x-y-diagram, it doesn't matter if the values are scaled in any way.
        # After all, the tranlation to XYZ is linear, and then to xyY it's (X/(X+Y+Z),
        # Y/(X+Y+Z), Y), so the factor will only be present in the last component which
        # is discarded. To make the plot a bit brighter, scale the colors up as much as
        # possible.
        rgb_linear /= numpy.max(rgb_linear, axis=0)

    srgb_linear = SrgbLinear()
    xyz = srgb_linear.to_xyz100(rgb_linear)
    xyy_vals = xy_to_2d(_xyy_from_xyz100(xyz)[:2])

    # Unfortunately, one cannot use tripcolors with explicit RGB specification
    # (see <https://github.com/matplotlib/matplotlib/issues/10265>). As a
    # workaround, associate range(n) data with the points and create a colormap
    # that associates the integer values with the respective RGBs.
    z = numpy.arange(xyy_vals.shape[1])
    rgb = srgb_linear.to_rgb1(rgb_linear)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "gamut", rgb.T, N=len(rgb.T)
    )

    triang = matplotlib.tri.Triangulation(xyy_vals[0], xyy_vals[1])
    plt.tripcolor(triang, z, shading="gouraud", cmap=cmap)
