import colorio


def test_rgb_gamut():
    colorspace = colorio.cs.CIELAB()
    colorio.plot_rgb_gamut(colorspace, n=51, show_grid=True)
    # p.show()
