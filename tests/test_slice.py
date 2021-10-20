import colorio


def test_rgb_slice():
    colorio.plot_rgb_slice(
        colorio.cs.CIELAB(), lightness=50.0, camera_elevation=10.0, off_screen=True
    )
