"""
https://en.wikipedia.org/wiki/CIECAM02
"""
import npx
import numpy as np
from numpy.typing import ArrayLike

from .._exceptions import ColorioError
from ..cat import cat02
from ..illuminants import whitepoints_cie1931
from ._color_space import ColorSpace
from ._helpers import register

M_hpe = np.array(
    [
        [+0.38971, +0.68898, -0.07868],
        [-0.22981, +1.18340, +0.04641],
        [+0.00000, +0.00000, +1.00000],
    ]
)


def compute_from(rgb_, cs):
    # Step 4: Calculate the post-adaptation cone response (resulting in dynamic range
    #         compression)
    alpha = (
        np.full(rgb_.shape, np.inf)
        if cs.F_L == np.inf
        else (cs.F_L * abs(rgb_) / 100) ** 0.42
    )
    # Omit the 0.1 here; that's canceled out in almost all cases below anyways (except
    # the computation of `t`).

    # Deal with alpha == 0, alpha == inf
    beta = np.empty(alpha.shape)
    idx0 = alpha < 1.0
    beta[idx0] = alpha[idx0] / (alpha[idx0] + 27.13)  # + 0.1
    idx1 = ~idx0
    beta[idx1] = 1.0 / (1.0 + 27.13 / alpha[idx1])  # + 0.1
    rgb_a_ = np.sign(rgb_) * 400 * beta  # + 0.1

    # Mix steps 5, 7, and part of step 10 here in one big dot-product.
    # Step 5: Calculate Redness-Greenness (a) , Yellowness-Blueness (b)
    #         components and hue angle (h)
    # Step 7: Calculate achromatic response A
    a, b, p2_, u = npx.dot(
        np.array([[11, -12, 1], [1, 1, -2], [40, 20, 1], [20, 20, 21]]), rgb_a_
    )
    a /= 11
    b /= 9
    p2_ /= 20
    u /= 20

    A = p2_ * cs.N_bb
    if np.any(A < 0):
        raise ColorioError("CIECAM02 breakdown")

    # Make sure that h is in [0, 360]
    h = np.rad2deg(np.arctan2(b, a)) % 360

    # Step 6: Calculate eccentricity (e_t) and hue composition (H), using the unique hue
    #         data given in Table 2.4.
    h_ = (h - cs.h[0]) % 360 + cs.h[0]
    e_t = (np.cos(np.deg2rad(h_) + 2) + 3.8) / 4
    # For searchsorted, we don't need the last entry of cs.h. That's because the method
    # will return len(arr)+1 if it wasn't successful. This also prevents difficulties
    # with nans.
    i = np.searchsorted(cs.h[:-1], h_) - 1
    beta = (h_ - cs.h[i]) * cs.e[i + 1]
    H = cs.H[i] + 100 * beta / (beta + cs.e[i] * (cs.h[i + 1] - h_))

    # Step 8: Calculate the correlate of lightness
    J = 100 * (A / cs.A_w) ** (cs.c * cs.z)

    # Step 9: Calculate the correlate of brightness
    sqrt_J_100 = np.sqrt(J / 100)
    Q = (4 / cs.c) * sqrt_J_100 * (cs.A_w + 4) * cs.F_L**0.25

    # Step 10: Calculate the correlates of chroma (C), colourfulness (M) and saturation
    # (s)
    #
    # Note the extra 0.305 here from the adaptation in rgb_a_ above.
    p1_ = 50000 / 13 * e_t * cs.N_c * cs.N_cb
    t = p1_ * np.hypot(a, b) / (u + 0.305)

    if np.any(t < 0):
        raise ColorioError("CIECAM02 breakdown")

    alpha = t**0.9 * (1.64 - 0.29**cs.n) ** 0.73
    C = alpha * sqrt_J_100

    M = np.zeros(C.shape) if cs.F_L == np.inf else C * cs.F_L**0.25

    # ENH avoid division by Q=0 here.
    # s = 100 * np.sqrt(M/Q)
    s = 50 * np.sqrt(cs.c * alpha / (cs.A_w + 4))

    return np.array([J, C, H, h, M, s, Q])


def compute_to(data, description, cs):
    if description[0] == "J":
        J = data[0]
        # Q perhaps needed for C
        # Q = (4 / cs.c) * np.sqrt(J / 100) * (cs.A_w + 4) * cs.F_L ** 0.25
    else:
        # Step 1-1: Compute J from Q (if start from Q)
        assert description[0] == "Q"
        Q = data[0]
        J = 6.25 * (cs.c * Q / (cs.A_w + 4) / cs.F_L**0.25) ** 2

    # Step 1-2: Calculate t from C, M, or s
    if description[1] in ["C", "M"]:
        if description[1] == "M":
            M = data[1]
            C = M / cs.F_L**0.25
        else:
            C = data[1]

        # If C or M is given and equal 0, the value of `t` cannot algebraically deduced
        # just by C or M. However, from other considerations we know that it must be 0.
        alpha = np.zeros_like(J)
        np.divide(C, np.sqrt(J / 100), out=alpha, where=J != 0.0)
    else:
        assert description[1] == "s"
        s = data[1] / 100
        # C = s * s * Q / cs.F_L ** 0.25
        alpha = 4 * s * s * (cs.A_w + 4) / cs.c

    t = (alpha / (1.64 - 0.29**cs.n) ** 0.73) ** (1 / 0.9)

    if description[2] == "h":
        h = data[2]
    else:
        assert description[2] == "H"
        # Step 1-3: Calculate h from H (if start from H)
        H = data[2]
        i = np.searchsorted(cs.H, H) - 1
        Hi = cs.H[i]
        hi, hi1 = cs.h[i], cs.h[i + 1]
        ei, ei1 = cs.e[i], cs.e[i + 1]
        h_ = ((H - Hi) * (ei1 * hi - ei * hi1) - 100 * hi * ei1) / (
            (H - Hi) * (ei1 - ei) - 100 * ei1
        )
        h = np.mod(h_, 360)

    # Step 2: Calculate t, et , p1, p2 and p3
    e_t = 0.25 * (np.cos(h * np.pi / 180 + 2) + 3.8)
    A = cs.A_w * (J / 100) ** (1 / cs.c / cs.z)

    # no 0.305
    p2_ = A / cs.N_bb

    # Step 3: Calculate a and b
    # ENH Much more straightforward computation of a, b
    p1_ = e_t * 50000 / 13 * cs.N_c * cs.N_cb
    sinh = np.sin(h * np.pi / 180)
    cosh = np.cos(h * np.pi / 180)
    a, b = np.array([cosh, sinh]) * (
        23 * (p2_ + 0.305) * t / (23 * p1_ + 11 * t * cosh + 108 * t * sinh)
    )

    # Step 4: Calculate RGB_a_
    rgb_a_ = (
        npx.dot(
            np.array([[460, 451, 288], [460, -891, -261], [460, -220, -6300]]),
            np.array([p2_, a, b]),
        )
        / 1403
    )

    # Step 5: Calculate RGB_
    t = np.array(
        [
            1.0
            if cs.F_L == np.inf
            else ((27.13 * abs(r)) / (400 - abs(r))) ** (1 / 0.42) / cs.F_L
            for r in rgb_a_
        ]
    )
    rgb_ = np.sign(rgb_a_) * 100 * t

    return rgb_


class CIECAM02:
    """
    Ming Ronnier Luo and Changjun Li,
    CIECAM02 and Its Recent Developments,
    Chapter 2
    <https://link.springer.com/chapter/10.1007%2F978-1-4419-6190-7_2>
    <https://www.springer.com/cda/content/document/cda_downloaddocument/9781441961891-c1.pdf?SGWID=0-0-45-1337202-p173968189>

    Appendix: CIE Colour Appearance Model: CIECAM02
    Part 1: The Forward Mode

    c ..... surround parameter (average: 0.69, dim: 0.59, dark: 0.535)
    Y_b ... relative luminance of the background in %, i.e., 100 * Lb / Lw
            where Lb is the background luminance and Lw is the white luminance.
            There is some discussion on the value of Y_b, see
            <https://groups.google.com/g/sci.engr.color/c/3P9DXaCAWAI>
            <https://rawpedia.rawtherapee.com/CIECAM02/de>
            The suggestion is to use 18 or 20 (gray world theory).
    L_A ... luminance of the adapting field in cd/m^2.
            Can be approximated by Lw * Y_b / 100, or Lb.

    From the above article:
    Table 2.1 Parameter settings for some typical applications:

        * Surface color evaluation in a light booth:
           - Ambient illumination in lux (cd/m2): 1000 (318.3)
           - Scene or device white luminance: 318.3 cd/m2
           - L_A: 60 cd/m2
           - Adopted white point: Light booth
           - Surround ratio: 1
           - Surround: Average (c = 0.69)

        * Viewing self-luminous display at home:
           - Ambient illumination in lux (cd/m2): 38 (12)
           - Scene or device white luminance: 80 cd/m2
           - L_A: 20 cd/m2
           - Adopted white point: Display and ambient
           - Surround ratio: 0.15
           - Surround: Dim (c = 0.59)

        * Viewing slides in dark room:
           - Ambient illumination in lux (cd/m2): 0 (0)
           - Scene or device white luminance: 150 cd/m2
           - L_A: 30 cd/m2
           - Adopted white point: Projector
           - Surround ratio: 0
           - Surround: Dark (c = 0.525)

        * Viewing self-luminous display under office illumination:
           - Ambient illumination in lux (cd/m2): 500 (159.2)
           - Scene or device white luminance: 80 cd/m2
           - L_A: 15 cd/m2
           - Adopted white point: Display
           - Surround ratio: 2
           - Surround: Average (c = 0.69)

    Publication CIE 159:
    A colour appearance model for colour management systems: CIECAM02,
    <DOI: 10.1002/col.20198>.
    """

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

        # Nc and F are modelled as a function of c, and can be linearly interpolated.
        #
        # Surround
        # condition  Surround ratio  F       c       Nc      Application
        # -----------------------------------------------------------------------
        # Average    SR > 0.15       1.0     0.69    1.0     Viewing surface colors
        # Dim        0 < SR < 0.15   0.9     0.59    0.9     Viewing television
        # Dark       SR = 0          0.8     0.525   0.8     Projector in a dark room
        #
        # SR = Lsw / Ldw:
        #     ratio of the absolute luminance of the reference white (white point)
        #     measured in the surround field to the display area. The 0.2 coefficient
        #     derives from the "gray world" assumption (~18%–20% reflectivity). It tests
        #     whether the surround luminance is darker or brighter than medium gray.
        # F:  factor determining degree of adaptation
        # c:  impact of surrounding
        # Nc: chromatic induction factor

        c_vals = [0.525, 0.59, 0.69]
        F_Nc_vals = [0.8, 0.9, 1.0]
        assert 0.525 <= c <= 0.69, f"c = {c}"
        F = np.interp(c, c_vals, F_Nc_vals)
        self.c = c
        self.N_c = F

        self.M, self.Minv = cat02(
            whitepoint_source=whitepoint,
            whitepoint_target=[100.0, 100.0, 100.0],
            F=F,
            L_A=L_A,
        )
        self.M = M_hpe @ self.M
        self.Minv = self.Minv @ np.linalg.inv(M_hpe)

        k = 1 / (5 * L_A + 1)
        k4 = k**4
        l4 = 1 - k4
        self.F_L = k4 * L_A + 0.1 * l4**2 * np.cbrt(5 * L_A)

        self.n = Y_b / Y_w
        self.z = 1.48 + np.sqrt(self.n)
        self.N_bb = 0.725 / self.n**0.2
        self.N_cb = self.N_bb

        RGB_w_ = self.M @ whitepoint

        alpha = (self.F_L * RGB_w_ / 100) ** 0.42
        RGB_aw_ = 400 * alpha / (alpha + 27.13)
        self.A_w = np.dot([2, 1, 1 / 20], RGB_aw_) * self.N_bb

        self.h = np.array([20.14, 90.00, 164.25, 237.53, 380.14])
        self.e = np.array([0.8, 0.7, 1.0, 1.2, 0.8])
        self.H = np.array([0.0, 100.0, 200.0, 300.0, 400.0])

    def from_xyz100(self, xyz):
        # Step 1: Calculate (sharpened) cone responses (transfer
        #         colour-matching functions to sharper sensors)
        #
        # Step 2: Calculate the corresponding (sharpened) cone response
        #         (considering various luminance level and surround conditions
        #         included in D; hence, in DR, DG and DB)
        #
        # Step 3: Calculate the Hunt-Pointer-Estevez response

        # cat02: Illuminant color adaptation

        rgb_ = npx.dot(self.M, xyz)

        # Steps 4-10
        return compute_from(rgb_, self)

    def to_xyz100(self, data, description):
        """Input: J or Q; C, M or s; H or h"""
        # Steps 1-5
        rgb_ = compute_to(data, description, self)

        # Step 6: Calculate RC, GC and BC
        # rgb_c = dot(M_cat02, solve(M_hpe, rgb_))
        #
        # Step 7: Calculate R, G and B
        # rgb = (rgb_c.T / self.D_RGB).T
        #
        # Step 8: Calculate X, Y and Z
        # xyz = solve(M_cat02, rgb)
        return npx.dot(self.Minv, rgb_)


class CAM02(ColorSpace):
    labels = ("J'", "a'", "b'")
    k0 = 0

    def __init__(
        self,
        variant: str,
        c: float,
        Y_b: float,
        L_A: float,
        whitepoint: ArrayLike = whitepoints_cie1931["D65"],
    ):
        params = {
            "LCD": (0.77, 0.007, 0.0053),
            "SCD": (1.24, 0.007, 0.0363),
            "UCS": (1.00, 0.007, 0.0228),
        }
        self.K_L, self.c1, self.c2 = params[variant]
        self.ciecam02 = CIECAM02(c, Y_b, L_A, whitepoint)
        self.name = f"CAM02 ({variant})"

    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        J, _, _, h, M, _, _ = self.ciecam02.from_xyz100(xyz)
        J_ = (1 + 100 * self.c1) * J / (1 + self.c1 * J)
        M_ = 1 / self.c2 * np.log(1 + self.c2 * M)
        h_ = h / 180 * np.pi
        return np.array([J_, M_ * np.cos(h_), M_ * np.sin(h_)])

    def to_xyz100(self, jab: ArrayLike) -> np.ndarray:
        J_, a, b = np.asarray(jab)
        J = J_ / (1 - (J_ - 100) * self.c1)
        h = np.mod(np.arctan2(b, a), 2 * np.pi) / np.pi * 180
        M_ = np.hypot(a, b)
        M = (np.exp(M_ * self.c2) - 1) / self.c2
        return self.ciecam02.to_xyz100(np.array([J, M, h]), "JMh")


class CAM02LCD(CAM02):
    name = "CAM02 (LCD)"

    def __init__(self, c, Y_b, L_A, whitepoint):
        super().__init__("LCD", c, Y_b, L_A, whitepoint)


class CAM02SCD(CAM02):
    name = "CAM02 (SCD)"

    def __init__(self, c, Y_b, L_A, whitepoint):
        super().__init__("SCD", c, Y_b, L_A, whitepoint)


class CAM02UCS(CAM02):
    name = "CAM02 (UCS)"

    def __init__(self, c, Y_b, L_A, whitepoint):
        super().__init__("UCS", c, Y_b, L_A, whitepoint)


register("cam02lcd", CAM02LCD(0.69, 18.0, 20.0, whitepoints_cie1931["D65"]))
register("cam02scd", CAM02SCD(0.69, 18.0, 20.0, whitepoints_cie1931["D65"]))
register("cam02ucs", CAM02UCS(0.69, 18.0, 20.0, whitepoints_cie1931["D65"]))
