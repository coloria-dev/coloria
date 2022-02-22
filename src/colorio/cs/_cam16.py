import npx
import numpy as np
from numpy.typing import ArrayLike

from ..cat import cat16
from ..illuminants import whitepoints_cie1931
from ._ciecam02 import compute_from, compute_to
from ._color_space import ColorSpace
from ._helpers import register


class CAM16:
    """
    Li C, Li Z, Wang Z, et al.,
    Comprehensive color solutions: CAM16, CAT16, and CAM16-UCS.
    Color Res Appl. 2017;00:1-12.
    <https://doi.org/10.1002/col.22131>.

    For an explanation of the parameter, see CIECAM02.
    """

    name = "CAM16"

    def __init__(
        self,
        c: float,
        Y_b: float,
        L_A: float,
        whitepoint: ArrayLike = whitepoints_cie1931["D65"],
    ):
        # step0: Calculate all values/parameters which are independent of input
        #        samples
        whitepoint = np.asarray(whitepoint)
        Y_w = whitepoint[1]

        # Nc and F are modelled as a function of c, and can be linearly
        # interpolated.
        c_vals = [0.525, 0.59, 0.69]  # 0.525 vs. 0.535 in CIECAM02
        F_Nc_vals = [0.8, 0.9, 1.0]
        assert 0.525 <= c <= 0.69
        F = np.interp(c, c_vals, F_Nc_vals)
        self.c = c
        self.N_c = F

        self.M, self.Minv = cat16(
            whitepoint,
            whitepoint_target=[100.0, 100.0, 100.0],
            F=F,
            L_A=L_A,
            # Skip transformation back in to XYZ space because the the lightness
            # adaptation also happens in the transformed space. The fact that chromic
            # and lightness adaption can happen in the same space is one of the main
            # claims of CAM16.
            include_back_transform=False,
        )

        k = 1 / (5 * L_A + 1)
        l4 = 1 - k**4
        k4_L_A = 0.0 if L_A == np.inf else k**4 * L_A
        self.F_L = k4_L_A + 0.1 * l4**2 * np.cbrt(5 * L_A)

        self.n = Y_b / Y_w
        self.z = 1.48 + np.sqrt(self.n)
        self.N_bb = 0.725 / self.n**0.2
        self.N_cb = self.N_bb

        RGB_wc = self.M @ whitepoint

        alpha = (self.F_L * RGB_wc / 100) ** 0.42
        RGB_aw_ = np.array(
            [400.1 if a == np.inf else 400 * a / (a + 27.13) + 0.1 for a in alpha]
        )
        self.A_w = (np.dot([2, 1, 1 / 20], RGB_aw_) - 0.305) * self.N_bb

        self.h = np.array([20.14, 90.00, 164.25, 237.53, 380.14])
        self.e = np.array([0.8, 0.7, 1.0, 1.2, 0.8])
        self.H = np.array([0.0, 100.0, 200.0, 300.0, 400.0])

    def from_xyz100(self, xyz):
        # Step 1: Calculate 'cone' responses
        # rgb = dot(self.M16, xyz)
        # Step 2: Complete the color adaptation of the illuminant in
        #         the corresponding cone response space
        # rgb_c = (rgb.T * self.D_RGB).T
        rgb_c = npx.dot(self.M, xyz)
        return compute_from(rgb_c, self)

    def to_xyz100(self, data, description):
        """Input: J or Q; C, M or s; H or h"""
        rgb_c = compute_to(data, description, self)
        # Step 6: Calculate R, G and B
        # rgb = (rgb_c.T / self.D_RGB).T
        # Step 7: Calculate X, Y and Z
        # xyz = self.solve_M16(rgb)
        return npx.dot(self.Minv, rgb_c)


class CAM16UCS(ColorSpace):
    name = "CAM16 (UCS)"
    labels = ("J'", "a'", "b'")
    k0 = 0

    def __init__(
        self,
        c: float,
        Y_b: float,
        L_A: float,
        whitepoint: ArrayLike = whitepoints_cie1931["D65"],
    ):
        self.K_L = 1.0
        self.c1 = 0.007
        self.c2 = 0.0228
        self.cam16 = CAM16(c, Y_b, L_A, whitepoint)

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        J, _, _, h, M, _, _ = self.cam16.from_xyz100(xyz)
        J_ = (1 + 100 * self.c1) * J / (1 + self.c1 * J)
        M_ = 1 / self.c2 * np.log(1 + self.c2 * M)
        h_ = h * np.pi / 180
        return np.array([J_, M_ * np.cos(h_), M_ * np.sin(h_)])

    def to_xyz100(self, jab: ArrayLike) -> np.ndarray:
        J_, a, b = np.asarray(jab)
        J = J_ / (1 - (J_ - 100) * self.c1)
        h = np.mod(np.arctan2(b, a) / np.pi * 180, 360)
        M_ = np.hypot(a, b)
        M = (np.exp(M_ * self.c2) - 1) / self.c2
        return self.cam16.to_xyz100(np.array([J, M, h]), "JMh")


register("cam16ucs", CAM16UCS(0.69, 18.0, 20.0))
