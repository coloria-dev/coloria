import numpy

from ._color_space import ColorSpace
from ._linalg import dot
from .illuminants import whitepoints_cie1931


class OsaUcs(ColorSpace):
    """
    David L. MacAdam,
    Uniform color scales,
    Journal of the Optical Society of America,
    Volume 64, Number 12, December 1974,
    <https://doi.org/10.1364/JOSAA.24.001823>,
    <https://en.wikipedia.org/wiki/OSA-UCS>.
    """

    def __init__(self, whitepoint=whitepoints_cie1931["D65"]):
        self.labels = ["L", "g", "j"]

        self.M = numpy.array(
            [
                [+0.7990, 0.4194, -0.1648],
                [-0.4493, 1.3265, +0.0927],
                [-0.1149, 0.3394, +0.7170],
            ]
        )
        self.Minv = numpy.linalg.inv(self.M)
        return

    def from_xyz100(self, xyz100):
        X, Y, Z = xyz100
        sum_xyz = numpy.sum(xyz100, axis=0)
        x = X / sum_xyz
        y = Y / sum_xyz

        K = (
            4.4934 * x ** 2
            + 4.3034 * y ** 2
            - 4.276 * x * y
            - 1.3744 * x
            - 2.5643 * y
            + 1.8103
        )
        Y0 = Y * K

        #  L' is L in original article
        L_prime = 5.9 * (numpy.cbrt(Y0) - 2 / 3 + 0.042 * numpy.cbrt(Y0 - 30))

        C = L_prime / (5.9 * (numpy.cbrt(Y0) - 2 / 3))
        R, G, B = dot(self.M, xyz100)

        a = -13.7 * numpy.cbrt(R) + 17.7 * numpy.cbrt(G) - 4 * numpy.cbrt(B)
        b = 1.7 * numpy.cbrt(R) + 8 * numpy.cbrt(G) - 9.7 * numpy.cbrt(B)

        L = (L_prime - 14.3993) / numpy.sqrt(2)
        g = C * a
        j = C * b

        return numpy.array([L, g, j])

    def to_xyz100(self, lgj, tol=1.0e-13, max_num_newton_steps=100):
        # Renbo Cao, H. Joel Trussell, and Renzo Shamey,
        # Comparison of the performance of inverse transformation methods from OSA-UCS
        # to CIEXYZ,
        # J Opt Soc Am A Opt Image Sci Vis., 2013 Aug 1,30(8):1508-15,
        # <https://doi.org/10.1364/JOSAA.30.001508>.
        L, g, j = lgj

        L_prime = L * numpy.sqrt(2) + 14.3993

        # Use Newton to find cbrt(Y0) =: t
        # f(t) = (L' / 5.9 + 2/3 - t) ** 3 - 0.042**3 * (t**3 - 30)
        # df/dt = -3 * (L' / 5.9 + 2/3 - t) ** 2 - 0.042**3 * 3 * t**2
        #
        # The following might be a good initial guess since 0.042 ** 3 * (t ** 3 - 30)
        # is small. An alternative is simply t=0.
        t = L_prime / 5.9 + 2 / 3
        ft = (L_prime / 5.9 + 2 / 3 - t) ** 3 - 0.042 ** 3 * (t ** 3 - 30)
        k = 0
        while numpy.any(numpy.abs(ft) > tol):
            if k >= max_num_newton_steps:
                raise RuntimeError(
                    "OSA-USC.to_xyz100 exceeded max number of Newton steps"
                )
            dfdt = -3 * (L_prime / 5.9 + 2 / 3 - t) ** 2 - 0.042 ** 3 * 3 * t ** 2
            t -= ft / dfdt
            ft = (L_prime / 5.9 + 2 / 3 - t) ** 3 - 0.042 ** 3 * (t ** 3 - 30)
            k += 1

        Y0 = t ** 3
        C = L_prime / (5.9 * (t - 2 / 3))
        a = g / C
        b = j / C

        # a = -13.7 * numpy.cbrt(R) + 17.7 * numpy.cbrt(G) - 4 * numpy.cbrt(B)
        # b = 1.7 * numpy.cbrt(R) + 8 * numpy.cbrt(G) - 9.7 * numpy.cbrt(B)
        # A = [[17.7, -4], [8, -9.7]]
        det = -17.7 * 9.7 + 4 * 8
        ap = (-9.7 * a + 4 * b) / det
        bp = (-8 * a + 17.7 * b) / det
        # inv = [[-9.7, 4], [-8, 17.7]] / det

        def f_df(omega):
            cbrt_RGB = numpy.array([omega, omega + ap, omega + bp])
            # print(cbrt_RGB)
            # cbrt_R = omega
            # ap = (a + 13.7 * omega) / det
            # bp = (b - 1.7 * omega) / det
            # cbrt_G = -9.7 * ap + 4 * bp
            # cbrt_B = -8 * ap + 17.7 * bp

            # A = numpy.array([[-13.7, 17.7, -4], [1.7, 8, -9.7], [0.0, 1.0, 0.0]])
            # Ainv = numpy.linalg.inv(A)
            # rhs = numpy.array(
            #     [numpy.full(omega.shape, a), numpy.full(omega.shape, b), omega]
            # )
            # cbrt_RGB = dot(Ainv, rhs)
            # # a = -13.7 * numpy.cbrt(R) + 17.7 * numpy.cbrt(G) - 4 * numpy.cbrt(B)
            # # b = 1.7 * numpy.cbrt(R) + 8 * numpy.cbrt(G) - 9.7 * numpy.cbrt(B)

            # print(cbrt_RGB)
            RGB = cbrt_RGB ** 3
            # print(RGB)
            xyz100 = dot(self.Minv, RGB)

            # evil = dot(self.Minv.T, [1, 1, 1])
            # print("evil", evil, numpy.sum(evil))
            # print(dot(evil, RGB))
            # print()

            # print(xyz100)

            X, Y, Z = xyz100
            sum_xyz = numpy.sum(xyz100, axis=0)
            # print(sum_xyz)
            x = X / sum_xyz
            y = Y / sum_xyz
            # print(x, y)
            K = (
                4.4934 * x ** 2
                + 4.3034 * y ** 2
                - 4.276 * x * y
                - 1.3744 * x
                - 2.5643 * y
                + 1.8103
            )
            # print(K)
            f = Y * K - Y0

            # df/domega
            # dcbrt_RGB = 1.0
            # dRGB = 3 * cbrt_RGB ** 2 * dcbrt_RGB
            dRGB = 3 * cbrt_RGB ** 2
            dxyz100 = dot(self.Minv, dRGB)

            dX, dY, dZ = dxyz100
            dsum_xyz = numpy.sum(dxyz100, axis=0)
            dx = (dX * sum_xyz - X * dsum_xyz) / sum_xyz ** 2
            dy = (dY * sum_xyz - Y * dsum_xyz) / sum_xyz ** 2
            dK = (
                4.4934 * 2 * x * dx
                + 4.3034 * 2 * y * dy
                - 4.276 * (dx * y + x * dy)
                - 1.3744 * dx
                - 2.5643 * dy
            )
            df = dY * K + Y * dK
            return f, df, xyz100

        # singular value of f: omega = 0.16865940093
        # this happens if x+y+z approx 0, some of them being negative
        # This means that the inverse function is not well-defined since it is not
        # bijective. :(
        #
        # print()
        # val, _, _ = f_df(numpy.array(1.747160437))
        # print("{:.6e}".format(val))
        # exit(1)
        #
        # x = numpy.linspace(-1.0, 20.0, 10000)
        # y, _, _ = f_df(x)
        # import matplotlib.pyplot as plt

        # plt.plot(x, y)
        # plt.grid()
        # # plt.ylim(-25, 500)
        # # plt.axes().set_aspect('equal')
        # plt.show()

        # KOBAYASI, Mituo* and YosIKI, Kayoko*
        # An Effective Conversion Algorithm from OSA-UCS to CIEXYZ,
        # one reads:
        #
        # > Examining the property of phi(w) for many (L,j,g)'s, it is found that the
        # > function (w) is monotone increasing, convex downward, and smooth.
        #
        # None of this is correct. Unfortunately, phi has a singularity and more than
        # one zeros. It is thus very important to have a reasonable initial guess.
        # Indeed, one reads further:
        #
        # > When Newton-Raphson method is applied to solve the equation (11) or (12), it
        # > is important to find a good estimation of an initial approximation of the
        # > solution. We utilize the analytical characteristics of the equation (e.g.
        # > smoothness, monotony, and convexity of a function) to find a good
        # > approximation for the exact solution, and to garantee the global convergence
        # > of the iteration. (Details are omitted.)
        #
        # Ah, the good old omittted details.
        # We find that it is crucial to start the Newton iteration to the right of the
        # solution; there, the function seems well-behaved. Hence, chose the value
        # cbrt_R max, with X=Y=100, Z=0.
        #
        omega = numpy.cbrt((0.7990 + 0.4194) * 100)
        fomega, dfdomega_val, xyz100 = f_df(omega)
        k = 0
        while numpy.any(numpy.abs(fomega) > tol):
            if k >= max_num_newton_steps:
                raise RuntimeError(
                    "OSA-USC.to_xyz100 exceeded max number of Newton steps"
                )
            omega -= fomega / dfdomega_val
            fomega, dfdomega_val, xyz100 = f_df(omega)
            k += 1

        return xyz100
