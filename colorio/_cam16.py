import numpy

from ._ciecam02 import compute_from, compute_to
from ._color_space import ColorSpace
from ._linalg import dot
from .illuminants import whitepoints_cie1931


class CAM16:
    """
    Li C, Li Z, Wang Z, et al.,
    Comprehensive color solutions: CAM16, CAT16, and CAM16-UCS.
    Color Res Appl. 2017;00:1-12.
    <https://doi.org/10.1002/col.22131>.
    """

    def __init__(
        self, c, Y_b, L_A, exact_inversion=True, whitepoint=whitepoints_cie1931["D65"]
    ):
        # step0: Calculate all values/parameters which are independent of input
        #        samples
        Y_w = whitepoint[1]

        # Nc and F are modelled as a function of c, and can be linearly
        # interpolated.
        c_vals = [0.525, 0.59, 0.69]  # 0.525 vs. 0.535 in CIECAM02
        F_Nc_vals = [0.8, 0.9, 1.0]
        assert 0.525 <= c <= 0.69
        F = numpy.interp(c, c_vals, F_Nc_vals)
        self.c = c
        self.N_c = F

        self.M16 = numpy.array(
            [
                [+0.401288, +0.650173, -0.051461],
                [-0.250268, +1.204414, +0.045854],
                [-0.002079, +0.048952, +0.953127],
            ]
        )
        RGB_w = numpy.dot(self.M16, whitepoint)

        D = F * (1 - 1 / 3.6 * numpy.exp((-L_A - 42) / 92))
        D = min(D, 1.0)
        D = max(D, 0.0)

        self.D_RGB = D * Y_w / RGB_w + 1 - D

        k = 1 / (5 * L_A + 1)
        l4 = 1 - k ** 4
        k4_L_A = 0.0 if L_A == numpy.inf else k ** 4 * L_A
        self.F_L = k4_L_A + 0.1 * l4 ** 2 * numpy.cbrt(5 * L_A)

        self.n = Y_b / Y_w
        self.z = 1.48 + numpy.sqrt(self.n)
        self.N_bb = 0.725 / self.n ** 0.2
        self.N_cb = self.N_bb

        RGB_wc = self.D_RGB * RGB_w
        alpha = (self.F_L * RGB_wc / 100) ** 0.42
        RGB_aw_ = numpy.array(
            [400.1 if a == numpy.inf else 400 * a / (a + 27.13) + 0.1 for a in alpha]
        )
        self.A_w = (numpy.dot([2, 1, 1 / 20], RGB_aw_) - 0.305) * self.N_bb

        self.h = numpy.array([20.14, 90.00, 164.25, 237.53, 380.14])
        self.e = numpy.array([0.8, 0.7, 1.0, 1.2, 0.8])
        self.H = numpy.array([0.0, 100.0, 200.0, 300.0, 400.0])

        self.M_ = (self.M16.T * self.D_RGB).T
        # The standard acutally recommends using this approximation as
        # inversion operation.
        approx_inv_M16 = numpy.array(
            [
                [+1.86206786, -1.01125463, +0.14918677],
                [+0.38752654, +0.62144744, -0.00897398],
                [-0.01584150, -0.03412294, +1.04996444],
            ]
        )
        self.invM_ = (
            numpy.linalg.inv(self.M_)
            if exact_inversion
            else approx_inv_M16 / self.D_RGB
        )
        return

    def from_xyz100(self, xyz):
        # Step 1: Calculate 'cone' responses
        # rgb = dot(self.M16, xyz)
        # Step 2: Complete the color adaptation of the illuminant in
        #         the corresponding cone response space
        # rgb_c = (rgb.T * self.D_RGB).T
        rgb_c = dot(self.M_, xyz)
        return compute_from(rgb_c, self)

    def to_xyz100(self, data, description):
        """Input: J or Q; C, M or s; H or h"""
        rgb_c = compute_to(data, description, self)
        # Step 6: Calculate R, G and B
        # rgb = (rgb_c.T / self.D_RGB).T
        # Step 7: Calculate X, Y and Z
        # xyz = self.solve_M16(rgb)
        return dot(self.invM_, rgb_c)


class CAM16UCS(ColorSpace):
    def __init__(
        self, c, Y_b, L_A, exact_inversion=True, whitepoint=whitepoints_cie1931["D65"]
    ):
        self.K_L = 1.0
        self.c1 = 0.007
        self.c2 = 0.0228
        self.cam16 = CAM16(c, Y_b, L_A, exact_inversion, whitepoint)
        self.labels = ["J'", "a'", "b'"]
        self.k0 = 0  # the index that corresponds to luminosity
        return

    def from_xyz100(self, xyz):
        J, C, H, h, M, s, Q = self.cam16.from_xyz100(xyz)
        J_ = (1 + 100 * self.c1) * J / (1 + self.c1 * J)
        M_ = 1 / self.c2 * numpy.log(1 + self.c2 * M)
        h_ = h * numpy.pi / 180
        return numpy.array([J_, M_ * numpy.cos(h_), M_ * numpy.sin(h_)])

    def to_xyz100(self, jab):
        J_, a, b = jab
        J = J_ / (1 - (J_ - 100) * self.c1)
        h = numpy.mod(numpy.arctan2(b, a) / numpy.pi * 180, 360)
        M_ = numpy.hypot(a, b)
        M = (numpy.exp(M_ * self.c2) - 1) / self.c2
        return self.cam16.to_xyz100(numpy.array([J, M, h]), "JMh")
