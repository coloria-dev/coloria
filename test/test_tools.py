import colorio


def test_flat_gamut():
    colorio.show_flat_gamut()


def test_xy_gamut_mesh():
    # points, cells =
    colorio.xy_gamut_mesh(0.05)
    # import meshio
    # meshio.write_points_cells("test.vtu", points, {"triangle": cells})


def test_srgb_gradient():
    # cs = colorio.cs.CIELAB()
    cs = colorio.cs.OKLAB()
    # cs = colorio.cs.DIN99()
    # cs = colorio.cs.CIELUV()
    cs.show_srgb_gradient([0, 0, 255], [255, 255, 0])


if __name__ == "__main__":
    test_srgb_gradient()
