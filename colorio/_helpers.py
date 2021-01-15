import numpy as np


def _find_Y(cs, xy, level, tol=1.0e-5):
    """Use bisection to find a matching Y value that projects the xy into the given
    level.
    """
    x, y = xy
    min_Y = 0.0
    xyz100 = np.array([min_Y / y * x, min_Y, min_Y / y * (1 - x - y)]) * 100
    min_val = cs.from_xyz100(xyz100)[cs.k0]
    assert min_val <= level

    # search for an appropriate max_Y to start with
    max_Y = 1.0
    while True:
        xyz100 = np.array([max_Y / y * x, max_Y, max_Y / y * (1 - x - y)]) * 100
        max_val = cs.from_xyz100(xyz100)[cs.k0]
        if max_val >= level:
            break
        max_Y *= 2

    while True:
        Y = (max_Y + min_Y) / 2
        xyz100 = np.array([Y / y * x, Y, Y / y * (1 - x - y)]) * 100
        val = cs.from_xyz100(xyz100)
        if abs(val[cs.k0] - level) < tol:
            break
        elif val[cs.k0] > level:
            max_Y = Y
        else:
            assert val[cs.k0] < level
            min_Y = Y

    return val
