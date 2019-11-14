import numpy

from ._color_space import ColorSpace
from ._linalg import dot, solve
from .illuminants import whitepoints_cie1931


class Osa(ColorSpace):
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

    def to_xyz100(self, lgj, tol=1.0e-13):
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
        # This might be a good initial guess since 0.042 ** 3 * (t ** 3 - 30) is small.
        # An alternative is simply t=0.
        t = L_prime / 5.9 + 2 / 3
        ft = (L_prime / 5.9 + 2 / 3 - t) ** 3 - 0.042 ** 3 * (t ** 3 - 30)
        while numpy.any(numpy.abs(ft) > tol):
            dfdt = -3 * (L_prime / 5.9 + 2 / 3 - t) ** 2 - 0.042 ** 3 * 3 * t ** 2
            t -= ft / dfdt
            ft = (L_prime / 5.9 + 2 / 3 - t) ** 3 - 0.042 ** 3 * (t ** 3 - 30)

        Y0 = t ** 3
        C = L_prime / (5.9 * (t - 2 / 3))
        a = g / C
        b = j / C

        # a = -13.7 * numpy.cbrt(R) + 17.7 * numpy.cbrt(G) - 4 * numpy.cbrt(B)
        # b = 1.7 * numpy.cbrt(R) + 8 * numpy.cbrt(G) - 9.7 * numpy.cbrt(B)
        # A = [[17.7, -4], [8, -9.7]]
        det = -17.7 * 9.7 + 4 * 8
        # inv = [[-9.7, 4], [-8, 17.7]] / det

        def f(omega):
            cbrt_R = omega
            ap = (a + 13.7 * omega) / det
            bp = (b - 1.7 * omega) / det
            cbrt_G = -9.7 * ap + 4 * bp
            cbrt_B = -8 * ap + 17.7 * bp
            # cbrt_G = (-9.7 * a + 4 * b) / det + (-9.7 * 13.7 - 4 * 1.7) / det * omega
            # cbrt_B = (-8 * a + 17.7 * b) / det + (-8 * 13.7 - 17.7 * 1.7) / det * omega

            RGB = numpy.array([cbrt_R, cbrt_G, cbrt_B]) ** 3
            xyz100 = dot(self.Minv, RGB)

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
            return Y * K - Y0

        def dfdomega(omega):
            cbrt_R = omega
            ap = (a + 13.7 * omega) / det
            bp = (b - 1.7 * omega) / det
            cbrt_G = -9.7 * ap + 4 * bp
            cbrt_B = -8 * ap + 17.7 * bp

            dcbrt_R = 1.0
            dap = 13.7 / det
            dbp = -1.7 / det
            dcbrt_G = -9.7 * dap + 4 * dbp
            dcbrt_B = -8 * dap + 17.7 * dbp

            RGB = numpy.array([cbrt_R, cbrt_G, cbrt_B]) ** 3
            dRGB = numpy.array([
                3 * cbrt_R**2 * dcbrt_R,
                3 * cbrt_G**2 * dcbrt_G,
                3 * cbrt_B**2 * dcbrt_B,
            ])
            xyz100 = dot(self.Minv, RGB)
            dxyz100 = dot(self.Minv, dRGB)

            X, Y, Z = xyz100
            dX, dY, dZ = dxyz100
            sum_xyz = numpy.sum(xyz100, axis=0)
            dsum_xyz = numpy.sum(dxyz100, axis=0)
            x = X / sum_xyz
            y = Y / sum_xyz
            dx = (dX * sum_xyz - X * dsum_xyz) / sum_xyz ** 2
            dy = (dY * sum_xyz - Y * dsum_xyz) / sum_xyz ** 2
            K = (
                4.4934 * x ** 2
                + 4.3034 * y ** 2
                - 4.276 * x * y
                - 1.3744 * x
                - 2.5643 * y
                + 1.8103
            )
            dK = (
                4.4934 * 2 * x * dx
                + 4.3034 * 2 * y * dy
                - 4.276 * (dx * y + x * dy)
                - 1.3744 * dx
                - 2.5643 * dy
            )
            return dY * K + Y * dK

        # singular value of f: omega = 0.16865940093
        # this happens if x+y+z approx 0, some of them being negative
        # This means that the inverse function is not well-defined since it is not
        # bijective. :(
        #
        # print()
        # print("{:.6e}".format(f(0.16865940093)))
        # exit(1)
        #
        # x = numpy.linspace(0.0, 5.0, 10000)
        # y = f(x)
        # import matplotlib.pyplot as plt
        # plt.plot(x, y)
        # plt.grid()
        # plt.show()
        # exit(1)

        # omega = 3.9981595815071427
        # omega = 3.99
        omega = 0.0
        fomega = f(omega)
        while numpy.any(numpy.abs(fomega) > tol):
            dfdomega_val = dfdomega(omega)
            omega -= fomega / dfdomega_val
            fomega = f(omega)

        # a + 13.7 * omega = + 17.7 * numpy.cbrt(G) - 4 * numpy.cbrt(B)
        # b - 1.7 * omega = 8 * numpy.cbrt(G) - 9.7 * numpy.cbrt(B)
        cbrt_R = omega
        ap = (a + 13.7 * omega) / det
        bp = (b - 1.7 * omega) / det
        cbrt_G = -9.7 * ap + 4 * bp
        cbrt_B = -8 * ap + 17.7 * bp

        RGB = numpy.array([cbrt_R, cbrt_G, cbrt_B]) ** 3
        xyz100 = dot(self.Minv, RGB)

        return xyz100
