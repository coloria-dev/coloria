import numpy as np

from colorio.cs import ColorCoordinates, convert


def test_conversion():
    color_coords = ColorCoordinates(
        np.array([[0.0, 1.0, 0.3], [1.5, 0.6, 0.5]]).T, "OKLAB"
    )

    # out = color_coords.convert("CIELAB")
    out = convert(color_coords, "CIELAB")

    print(out.data.T.tolist())

    ref = [
        [1.4816557743706664, 153.25609942581133, 258.2179692122982],
        [147.32553440561583, 235.88955497683463, 1447.9761970245452],
    ]
    assert np.all(np.abs(out.data.T - ref) < 1.0e-13 * (1.0 + np.abs(ref)))
