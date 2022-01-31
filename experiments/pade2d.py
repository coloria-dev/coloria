import numpy as np


def _create_tree(alpha, degree):
    return [
        np.array([alpha[d * (d + 1) // 2 + i] for i in range(d + 1)])
        for d in range(degree + 1)
    ]
    # return np.split(alpha, np.arange(1, degree+1))


def _get_xy_tree(xy, degree):
    """Evaluates the entire tree of 2d mononomials.

    The return value is a list of arrays, where `out[k]` hosts the `2*k+1`
    values of the `k`th level of the tree

        (0, 0)
        (1, 0)   (0, 1)
        (2, 0)   (1, 1)   (0, 2)
          ...      ...      ...
    """
    x, y = xy
    tree = [np.array([np.ones(x.shape, dtype=int)])]
    for d in range(degree):
        tree.append(np.concatenate([tree[-1] * x, [tree[-1][-1] * y]]))
    return tree


def _get_dx_tree(xy, degree):
    """
                                  0
                        1*(0, 0)  0
              2*(1, 0)  1*(0, 1)  0
    3*(2, 0)  2*(1, 1)  1*(0, 2)  0
      ...       ...      ...     ...
    """
    x, y = xy

    # build smaller tree
    one = np.array([np.ones(x.shape, dtype=int)])
    tree = [one]
    for d in range(1, degree):
        tree.append(
            np.concatenate(
                [
                    # Integer division `//` would be nice here, but
                    # <https://github.com/sympy/sympy/issues/14542>.
                    [tree[-1][0] / d * (d + 1) * x],
                    tree[-1] * y,
                ]
            )
        )

    # append zeros
    zero = np.array([np.zeros(x.shape, dtype=int)])
    tree = [zero] + [np.concatenate([t, zero]) for t in tree]
    return tree


def _get_dy_tree(xy, degree):
    """
     0
     0  1*(0, 0)
     0  1*(1, 0)  2*(0, 1)
     0  1*(2, 0)  2*(1, 1)  3*(0, 2)
    ...   ...       ...       ...
    """
    x, y = xy

    one = np.array([np.ones(x.shape, dtype=int)])
    tree = [one]
    for d in range(1, degree):
        tree.append(
            np.concatenate(
                [
                    tree[-1] * x,
                    # Integer division `//` would be nice here, but
                    # <https://github.com/sympy/sympy/issues/14542>.
                    [tree[-1][-1] / d * (d + 1) * y],
                ]
            )
        )

    # prepend zeros
    zero = np.array([np.zeros(x.shape, dtype=int)])
    tree = [zero] + [np.concatenate([zero, t]) for t in tree]
    return tree


class Pade2d:
    """Pad'e polynomial in from R^2 to R^2, i.e., both components are Pad'e
    functions in x and y. The function is characterized by four sets of
    coefficients, two for the components and two for numerator and denominator
    each.
    """

    def __init__(self, xy, degrees, ax, bx, ay, by):
        self.degrees = degrees
        self.set_coefficients(ax, bx, ay, by)
        self.set_xy(xy)
        return

    def set_coefficients(self, ax, bx, ay, by):
        assert len(ax) == (self.degrees[0] + 1) * (self.degrees[0] + 2) // 2
        assert len(bx) == (self.degrees[1] + 1) * (self.degrees[1] + 2) // 2
        assert len(ay) == (self.degrees[2] + 1) * (self.degrees[2] + 2) // 2
        assert len(by) == (self.degrees[3] + 1) * (self.degrees[3] + 2) // 2

        self.ax = ax
        self.ay = ay
        self.bx = bx
        self.by = by
        return

    def set_xy(self, xy):
        self.xy = xy

        xy_tree = _get_xy_tree(xy, max(self.degrees))
        self.xy_list = np.array([item for branch in xy_tree for item in branch])

        dx_tree = _get_dx_tree(xy, max(self.degrees))
        self.dx_list = np.array([item for branch in dx_tree for item in branch])

        dy_tree = _get_dy_tree(xy, max(self.degrees))
        self.dy_list = np.array([item for branch in dy_tree for item in branch])
        return

    def eval(self, xy=None):
        if xy is not None:
            self.set_xy(xy)

        ux = np.dot(self.ax, self.xy_list[: len(self.ax)])
        vx = np.dot(self.bx, self.xy_list[: len(self.bx)])
        uy = np.dot(self.ay, self.xy_list[: len(self.ay)])
        vy = np.dot(self.by, self.xy_list[: len(self.by)])

        return np.array([ux / vx, uy / vy])

    def jac(self, xy=None):
        """Get the Jacobian at (x, y)."""
        if xy is not None:
            self.set_xy(xy)

        ux = np.dot(self.ax, self.xy_list[: len(self.ax)])
        vx = np.dot(self.bx, self.xy_list[: len(self.bx)])
        uy = np.dot(self.ay, self.xy_list[: len(self.ay)])
        vy = np.dot(self.by, self.xy_list[: len(self.by)])

        ux_dx = np.dot(self.ax, self.dx_list[: len(self.ax)])
        vx_dx = np.dot(self.bx, self.dx_list[: len(self.bx)])
        uy_dx = np.dot(self.ay, self.dx_list[: len(self.ay)])
        vy_dx = np.dot(self.by, self.dx_list[: len(self.by)])

        ux_dy = np.dot(self.ax, self.dy_list[: len(self.ax)])
        vx_dy = np.dot(self.bx, self.dy_list[: len(self.bx)])
        uy_dy = np.dot(self.ay, self.dy_list[: len(self.ay)])
        vy_dy = np.dot(self.by, self.dy_list[: len(self.by)])

        jac = np.array(
            [
                [
                    (ux_dx * vx - vx_dx * ux) / vx**2,
                    (ux_dy * vx - vx_dy * ux) / vx**2,
                ],
                [
                    (uy_dx * vy - vy_dx * uy) / vy**2,
                    (uy_dy * vy - vy_dy * uy) / vy**2,
                ],
            ]
        )
        return jac

    def print(self):
        tree_ax = _create_tree(self.ax, self.degrees[0])
        tree_bx = _create_tree(self.bx, self.degrees[1])
        tree_ay = _create_tree(self.ay, self.degrees[2])
        tree_by = _create_tree(self.by, self.degrees[3])
        for k, vals in enumerate([tree_ax, tree_bx, tree_ay, tree_by]):
            for val in vals:
                print(val)
            print()
        return
