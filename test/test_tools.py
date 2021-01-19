import colorio


def test_flat_gamut():
    colorio.show_flat_gamut()


def test_xy_gamut_mesh():
    # points, cells =
    colorio.xy_gamut_mesh(0.05)
    # import meshio
    # meshio.write_points_cells("test.vtu", points, {"triangle": cells})
