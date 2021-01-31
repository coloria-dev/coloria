def _get_visible_gamut_mesh(colorspace, observer, illuminant):
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

    if not colorspace.is_origin_well_defined:
        values = values[1:]
        cells = cells[~np.any(cells == 0, axis=1)]
        cells -= 1

    pts = colorspace.from_xyz100(values.T).T
    return pts, cells


def save_visible_gamut(colorspace, observer, illuminant, filename):
    import meshio

    pts, cells = _get_visible_gamut_mesh(colorspace, observer, illuminant)
    meshio.write_points_cells(filename, pts, cells={"triangle": cells})


def show_visible_slice(*args, **kwargs):
    plt.figure()
    plot_visible_slice(*args, **kwargs)
    plt.show()
    plt.close()


def save_visible_slice(filename, *args, **kwargs):
    plt.figure()
    plot_visible_slice(*args, **kwargs)
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


def plot_visible_slice(colorspace, lightness, outline_prec=1.0e-2, fill_color="0.8"):
    # first plot the monochromatic outline
    mono_xy, conn_xy = get_mono_outline_xy(
        observer=observers.cie_1931_2(), max_stepsize=outline_prec
    )

    mono_vals = np.array([_find_Y(colorspace, xy, lightness) for xy in mono_xy])
    conn_vals = np.array([_find_Y(colorspace, xy, lightness) for xy in conn_xy])

    k1, k2 = [k for k in [0, 1, 2] if k != colorspace.k0]
    plt.plot(mono_vals[:, k1], mono_vals[:, k2], "-", color="k")
    plt.plot(conn_vals[:, k1], conn_vals[:, k2], ":", color="k")
    #
    if fill_color is not None:
        xyz = np.vstack([mono_vals, conn_vals[1:]])
        plt.fill(xyz[:, k1], xyz[:, k2], facecolor=fill_color, zorder=0)

    plt.axis("equal")
    plt.xlabel(colorspace.labels[k1])
    plt.ylabel(colorspace.labels[k2])
    plt.title(f"{colorspace.labels[colorspace.k0]} = {lightness}")
