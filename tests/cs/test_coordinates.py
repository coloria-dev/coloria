from colorio.cs import ColorCoordinates


def test_operations():
    coords = ColorCoordinates([1.0, 2.0, 3.0], "CIELAB")

    coords * 2
    2 * coords
    coords + 1
    1 + coords
    coords + coords

    coords < 0
    coords <= 0
    coords > 0
    coords >= 0
    coords == 0
