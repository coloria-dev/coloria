import npx
import numpy as np
from numpy.typing import ArrayLike

from ._color_space import ColorSpace
from ._helpers import register


class OsaUcs(ColorSpace):
    """
    David L. MacAdam,
    Uniform color scales,
    Journal of the Optical Society of America,
    Volume 64, Number 12, December 1974,
    <https://doi.org/10.1364/JOSA.64.001691>,
    <https://en.wikipedia.org/wiki/OSA-UCS>.

    Nico Schlömer,
    On the conversion from OSA-UCS to CIEXYZ,
    <https://arxiv.org/abs/1911.08323>.
    """

    name = "OSA-UCS"
    labels = ("L", "j", "g")
    k0 = 0

    def __init__(self):
        # The whitepoint used in the experiments (MacAdam 1974 color distances) was the
        # 10-degree observer under the D65 illuminant.

        self.M = np.array(
            [
                [+0.7990, 0.4194, -0.1648],
                [-0.4493, 1.3265, +0.0927],
                [-0.1149, 0.3394, +0.7170],
            ]
        )
        self.Minv = np.linalg.inv(self.M)

    def from_xyz100(self, xyz100: ArrayLike) -> np.ndarray:
        X, Y, _ = np.asarray(xyz100)
        s = np.sum(xyz100, axis=0)

        # Avoid division by s, could be 0.
        YKs2 = (
            4.4934 * Y * X**2
            + 4.3034 * Y**3
            - 4.276 * X * Y**2
            - 1.3744 * X * Y * s
            - 2.5643 * Y**2 * s
            + 1.8103 * Y * s**2
        )
        Y0 = np.zeros_like(s)
        np.divide(YKs2, s**2, out=Y0, where=s != 0.0)

        #  L' is L in original article
        L_prime = 5.9 * (np.cbrt(Y0) - 2 / 3 + 0.042 * np.cbrt(Y0 - 30))

        C = L_prime / (5.9 * (np.cbrt(Y0) - 2 / 3))
        R, G, B = npx.dot(self.M, xyz100)

        a = -13.7 * np.cbrt(R) + 17.7 * np.cbrt(G) - 4 * np.cbrt(B)
        b = 1.7 * np.cbrt(R) + 8 * np.cbrt(G) - 9.7 * np.cbrt(B)

        L = (L_prime - 14.3993) / np.sqrt(2)
        g = C * a
        j = C * b

        return np.array([L, j, g])

    def to_xyz100(
        self, ljg: ArrayLike, tol: float = 1.0e-13, max_num_newton_steps: int = 100
    ):
        # Renbo Cao, H. Joel Trussell, and Renzo Shamey,
        # Comparison of the performance of inverse transformation methods from OSA-UCS
        # to CIEXYZ,
        # J Opt Soc Am A Opt Image Sci Vis., 2013 Aug 1,30(8):1508-15,
        # <https://doi.org/10.1364/JOSAA.30.001508>.
        L, j, g = np.asarray(ljg)

        L_prime = L * np.sqrt(2) + 14.3993

        # Use Cardano to find cbrt(Y0) =: t
        # 0 = (L' / 5.9 + 2/3 - t) ** 3 - 0.042**3 * (t**3 - 30)
        # Note that the above function is monotonically decreasing and of third order,
        # hence has exactly one root.
        #
        u = L_prime / 5.9 + 2 / 3
        v = 0.042**3
        # Polynomial coefficients
        a = -(v + 1)
        b = 3 * u
        c = -3 * u**2
        d = u**3 + v * 30
        # val = a * t ** 3 + b * t ** 2 + c * t + d
        #
        # x = t + b / (3 * a)
        p = (3 * a * c - b**2) / (3 * a**2)
        q = (2 * b**3 - 9 * a * b * c + 27 * a**2 * d) / (27 * a**3)
        # val = (x ** 3 + p * x + q) * a
        #
        # No need to assert this: We already know from the original expression that the
        # equation has exactly one solution, so this must be >0.
        # assert np.all((p / 3) ** 3 + (q / 2) ** 2 > 0)
        #
        s = np.sqrt((q / 2) ** 2 + (p / 3) ** 3)
        t = np.cbrt(-q / 2 + s) + np.cbrt(-q / 2 - s)
        t -= b / (3 * a)

        Y0 = t**3
        C = L_prime / (5.9 * (t - 2 / 3))
        a = g / C
        b = j / C

        # a = -13.7 * np.cbrt(R) + 17.7 * np.cbrt(G) - 4 * np.cbrt(B)
        # b = 1.7 * np.cbrt(R) + 8 * np.cbrt(G) - 9.7 * np.cbrt(B)
        # A = [[17.7, -4], [8, -9.7]]
        det = -17.7 * 9.7 + 4 * 8
        ap = (-9.7 * a + 4 * b) / det
        bp = (-8 * a + 17.7 * b) / det
        # inv = [[-9.7, 4], [-8, 17.7]] / det

        def f_df(omega):
            omega = np.full_like(ap, omega)
            cbrt_RGB = np.array([omega, omega + ap, omega + bp])
            # A = np.array([[-13.7, 17.7, -4], [1.7, 8, -9.7], [0.0, 1.0, 0.0]])
            # Ainv = np.linalg.inv(A)
            # rhs = np.array(
            #     [np.full(omega.shape, a), np.full(omega.shape, b), omega]
            # )
            # cbrt_RGB = dot(Ainv, rhs)
            # # a = -13.7 * np.cbrt(R) + 17.7 * np.cbrt(G) - 4 * np.cbrt(B)
            # # b = 1.7 * np.cbrt(R) + 8 * np.cbrt(G) - 9.7 * np.cbrt(B)

            RGB = cbrt_RGB**3
            xyz100 = npx.dot(self.Minv, RGB)

            X, Y, _ = xyz100
            sum_xyz = np.sum(xyz100, axis=0)
            x = X / sum_xyz
            y = Y / sum_xyz
            K = (
                4.4934 * x**2
                + 4.3034 * y**2
                - 4.276 * x * y
                - 1.3744 * x
                - 2.5643 * y
                + 1.8103
            )
            f = Y * K - Y0

            # df/domega
            # dcbrt_RGB = 1.0
            # dRGB = 3 * cbrt_RGB ** 2 * dcbrt_RGB
            dRGB = 3 * cbrt_RGB**2
            dxyz100 = npx.dot(self.Minv, dRGB)

            dX, dY, _ = dxyz100
            dsum_xyz = np.sum(dxyz100, axis=0)
            dx = (dX * sum_xyz - X * dsum_xyz) / sum_xyz**2
            dy = (dY * sum_xyz - Y * dsum_xyz) / sum_xyz**2
            dK = (
                4.4934 * 2 * x * dx
                + 4.3034 * 2 * y * dy
                - 4.276 * (dx * y + x * dy)
                - 1.3744 * dx
                - 2.5643 * dy
            )
            df = dY * K + Y * dK
            return f, df, xyz100

        # In
        #   Kobayasi, Mituo and Yosiki, Kayoko
        #   An Effective Conversion Algorithm from OSA-UCS to CIEXYZ,
        # one reads:
        #
        # > Examining the property of psi(w) for many (L,j,g)'s, it is found that the
        # > function (w) is monotone increasing, convex downward, and smooth.
        #
        # None of this is correct. Unfortunately, psi has a singularity and more than
        # one zeros. It is thus very important to have a reasonable initial guess.
        # Indeed, one reads further:
        #
        # > When Newton-Raphson method is applied to solve the equation (11) or (12), it
        # > is important to find a good estimation of an initial approximation of the
        # > solution. We utilize the analytical characteristics of the equation (e.g.
        # > smoothness, monotony, and convexity of a function) to find a good
        # > approximation for the exact solution, and to guarantee the global
        # > convergence of the iteration. (Details are omitted.)
        #
        # Ah, the good old omittted details.
        # We find that it is crucial to start the Newton iteration to the right of the
        # solution; there, the function seems well-behaved. Hence, chose the value
        # cbrt_R max, with X=Y=100, Z=0.
        #
        omega = np.cbrt((0.7990 + 0.4194) * 100)
        fomega, dfdomega_val, xyz100 = f_df(omega)
        k = 0
        while np.any(np.abs(fomega) > tol):
            if k >= max_num_newton_steps:
                raise RuntimeError(
                    "OSA-USC.to_xyz100 exceeded max number of Newton steps"
                )
            omega -= fomega / dfdomega_val
            fomega, dfdomega_val, xyz100 = f_df(omega)
            k += 1

        return xyz100


register("osaucs", OsaUcs())
