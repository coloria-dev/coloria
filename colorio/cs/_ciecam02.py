import numpy

from .._exceptions import ColorioError
from .._linalg import dot
from ..illuminants import whitepoints_cie1931
from ._color_space import ColorSpace


def compute_from(rgb_, cs):
    # Step 4: Calculate the post-adaptation cone response (resulting in dynamic range
    #         compression)
    alpha = (
        numpy.full(rgb_.shape, numpy.inf)
        if cs.F_L == numpy.inf
        else (cs.F_L * abs(rgb_) / 100) ** 0.42
    )
    # Omit the 0.1 here; that's canceled out in almost all cases below anyways (except
    # the computation of `t`).

    # Deal with alpha == 0, alpha == inf
    beta = numpy.empty(alpha.shape)
    idx = alpha < 1.0
    beta[idx] = alpha[idx] / (alpha[idx] + 27.13)  # + 0.1
    idx = ~idx
    beta[idx] = 1.0 / (1.0 + 27.13 / alpha[idx])  # + 0.1
    rgb_a_ = numpy.sign(rgb_) * 400 * beta  # + 0.1

    # Mix steps 5, 7, and part of step 10 here in one big dot-product.
    # Step 5: Calculate Redness-Greenness (a) , Yellowness-Blueness (b)
    #         components and hue angle (h)
    # Step 7: Calculate achromatic response A
    a, b, p2_, u = dot(
        numpy.array([[11, -12, 1], [1, 1, -2], [40, 20, 1], [20, 20, 21]]), rgb_a_
    )
    a /= 11
    b /= 9
    p2_ /= 20
    u /= 20

    A = p2_ * cs.N_bb
    if numpy.any(A < 0):
        raise ColorioError("CIECAM02 breakdown")

    # Make sure that h is in [0, 360]
    h = numpy.rad2deg(numpy.arctan2(b, a)) % 360

    # Step 6: Calculate eccentricity (e_t) and hue composition (H), using the unique hue
    #         data given in Table 2.4.
    h_ = (h - cs.h[0]) % 360 + cs.h[0]
    e_t = (numpy.cos(numpy.deg2rad(h_) + 2) + 3.8) / 4
    i = numpy.searchsorted(cs.h, h_) - 1
    beta = (h_ - cs.h[i]) * cs.e[i + 1]
    H = cs.H[i] + 100 * beta / (beta + cs.e[i] * (cs.h[i + 1] - h_))

    # Step 8: Calculate the correlate of lightness
    J = 100 * (A / cs.A_w) ** (cs.c * cs.z)

    # Step 9: Calculate the correlate of brightness
    sqrt_J_100 = numpy.sqrt(J / 100)
    Q = (4 / cs.c) * sqrt_J_100 * (cs.A_w + 4) * cs.F_L ** 0.25

    # Step 10: Calculate the correlates of chroma (C), colourfulness (M) and saturation
    # (s)
    #
    # Note the extra 0.305 here from the adaptation in rgb_a_ above.
    p1_ = 50000 / 13 * e_t * cs.N_c * cs.N_cb
    t = p1_ * numpy.hypot(a, b) / (u + 0.305)

    alpha = t ** 0.9 * (1.64 - 0.29 ** cs.n) ** 0.73
    C = alpha * sqrt_J_100

    M = numpy.zeros(C.shape) if cs.F_L == numpy.inf else C * cs.F_L ** 0.25

    # ENH avoid division by Q=0 here.
    # s = 100 * numpy.sqrt(M/Q)
    s = 50 * numpy.sqrt(cs.c * alpha / (cs.A_w + 4))

    return numpy.array([J, C, H, h, M, s, Q])


def compute_to(data, description, cs):
    if description[0] == "J":
        J = data[0]
        # Q perhaps needed for C
        Q = (4 / cs.c) * numpy.sqrt(J / 100) * (cs.A_w + 4) * cs.F_L ** 0.25
    else:
        # Step 1-1: Compute J from Q (if start from Q)
        assert description[0] == "Q"
        Q = data[0]
        J = 6.25 * (cs.c * Q / (cs.A_w + 4) / cs.F_L ** 0.25) ** 2

    # Step 1-2: Calculate t from C, M, or s
    if description[1] in ["C", "M"]:
        if description[1] == "M":
            M = data[1]
            C = M / cs.F_L ** 0.25
        else:
            C = data[1]

        # If C or M is given and equal 0, the value of `t` cannot
        # algebraically deduced just by C or M. However, from other
        # considerations we know that it must be 0. Hence, allow division
        # by 0 and set nans to 0 afterwards.
        with numpy.errstate(invalid="ignore"):
            alpha = C / numpy.sqrt(J / 100)
        alpha = numpy.nan_to_num(alpha)
    else:
        assert description[1] == "s"
        s = data[1] / 100
        # C = s * s * Q / cs.F_L ** 0.25
        alpha = 4 * s * s * (cs.A_w + 4) / cs.c

    t = (alpha / (1.64 - 0.29 ** cs.n) ** 0.73) ** (1 / 0.9)

    if description[2] == "h":
        h = data[2]
    else:
        assert description[2] == "H"
        # Step 1-3: Calculate h from H (if start from H)
        H = data[2]
        i = numpy.searchsorted(cs.H, H) - 1
        Hi = cs.H[i]
        hi, hi1 = cs.h[i], cs.h[i + 1]
        ei, ei1 = cs.e[i], cs.e[i + 1]
        h_ = ((H - Hi) * (ei1 * hi - ei * hi1) - 100 * hi * ei1) / (
            (H - Hi) * (ei1 - ei) - 100 * ei1
        )
        h = numpy.mod(h_, 360)

    # Step 2: Calculate t, et , p1, p2 and p3
    e_t = 0.25 * (numpy.cos(h * numpy.pi / 180 + 2) + 3.8)
    A = cs.A_w * (J / 100) ** (1 / cs.c / cs.z)

    # no 0.305
    p2_ = A / cs.N_bb

    # Step 3: Calculate a and b
    # ENH Much more straightforward computation of a, b
    p1_ = e_t * 50000 / 13 * cs.N_c * cs.N_cb
    sinh = numpy.sin(h * numpy.pi / 180)
    cosh = numpy.cos(h * numpy.pi / 180)
    a, b = numpy.array([cosh, sinh]) * (
        23 * (p2_ + 0.305) * t / (23 * p1_ + 11 * t * cosh + 108 * t * sinh)
    )

    # Step 4: Calculate RGB_a_
    rgb_a_ = (
        dot(
            numpy.array([[460, 451, 288], [460, -891, -261], [460, -220, -6300]]),
            numpy.array([p2_, a, b]),
        )
        / 1403
    )

    # Step 5: Calculate RGB_
    t = numpy.array(
        [
            1.0
            if cs.F_L == numpy.inf
            else ((27.13 * abs(r)) / (400 - abs(r))) ** (1 / 0.42) / cs.F_L
            for r in rgb_a_
        ]
    )
    rgb_ = numpy.sign(rgb_a_) * 100 * t

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

    c ... surround parameter (average: 0.69, dim: 0.59, dark: 0.535)

    L_A ... in cd/m^2
    Table 2.1 Parameter settings for some typical applications:
    Surface colour evaluation in a light booth: 60
    Viewing self-luminous display at home: 20
    Viewing slides in dark room: 30
    Viewing self-luminous display under office illumination: 15

    Publication CIE 159:
    A colour appearance model for colour management systems: CIECAM02,
    <DOI: 10.1002/col.20198>.
    """

    def __init__(self, c, Y_b, L_A, whitepoint=whitepoints_cie1931["D65"]):
        # step0: Calculate all values/parameters which are independent of input
        #        samples
        Y_w = whitepoint[1]

        # Nc and F are modelled as a function of c, and can be linearly interpolated.
        c_vals = [0.525, 0.59, 0.69]
        F_Nc_vals = [0.8, 0.9, 1.0]
        assert 0.535 <= c <= 0.69
        F = numpy.interp(c, c_vals, F_Nc_vals)
        self.c = c
        self.N_c = F

        self.M_cat02 = numpy.array(
            [
                [+0.7328, +0.4296, -0.1624],
                [-0.7036, +1.6975, +0.0061],
                [+0.0030, +0.0136, +0.9834],
            ]
        )
        RGB_w = numpy.dot(self.M_cat02, whitepoint)

        D = F * (1 - 1 / 3.6 * numpy.exp((-L_A - 42) / 92))
        D = min(D, 1.0)
        D = max(D, 0.0)

        self.D_RGB = D * Y_w / RGB_w + 1 - D

        k = 1 / (5 * L_A + 1)
        k4 = k * k * k * k
        l4 = 1 - k4
        self.F_L = k4 * L_A + 0.1 * l4 * l4 * numpy.cbrt(5 * L_A)

        self.n = Y_b / Y_w
        self.z = 1.48 + numpy.sqrt(self.n)
        self.N_bb = 0.725 / self.n ** 0.2
        self.N_cb = self.N_bb

        RGB_wc = self.D_RGB * RGB_w

        self.M_hpe = numpy.array(
            [
                [+0.38971, +0.68898, -0.07868],
                [-0.22981, +1.18340, +0.04641],
                [+0.00000, +0.00000, +1.00000],
            ]
        )
        RGB_w_ = numpy.dot(self.M_hpe, numpy.linalg.solve(self.M_cat02, RGB_wc))

        alpha = (self.F_L * RGB_w_ / 100) ** 0.42
        RGB_aw_ = 400 * alpha / (alpha + 27.13)
        self.A_w = numpy.dot([2, 1, 1 / 20], RGB_aw_) * self.N_bb

        self.h = numpy.array([20.14, 90.00, 164.25, 237.53, 380.14])
        self.e = numpy.array([0.8, 0.7, 1.0, 1.2, 0.8])
        self.H = numpy.array([0.0, 100.0, 200.0, 300.0, 400.0])

        # Merge a bunch of matrices together here.
        self.M_ = numpy.dot(
            self.M_hpe,
            numpy.linalg.solve(self.M_cat02, (self.M_cat02.T * self.D_RGB).T),
        )
        # Alternative: LU decomposition. That introduces a scipy dependency
        # though and lusolve is slower than dot() as well.
        self.invM_ = numpy.linalg.inv(self.M_)

    def from_xyz100(self, xyz):
        # Step 1: Calculate (sharpened) cone responses (transfer
        #         colour-matching functions to sharper sensors)
        #
        # Step 2: Calculate the corresponding (sharpened) cone response
        #         (considering various luminance level and surround conditions
        #         included in D; hence, in DR, DG and DB)
        #
        # Step 3: Calculate the Hunt-Pointer-Estevez response
        rgb_ = dot(self.M_, xyz)
        # Steps 4-10
        return compute_from(rgb_, self)

    def to_xyz100(self, data, description):
        """Input: J or Q; C, M or s; H or h"""
        # Steps 1-5
        rgb_ = compute_to(data, description, self)

        # Step 6: Calculate RC, GC and BC
        # rgb_c = dot(self.M_cat02, solve(self.M_hpe, rgb_))
        #
        # Step 7: Calculate R, G and B
        # rgb = (rgb_c.T / self.D_RGB).T
        #
        # Step 8: Calculate X, Y and Z
        # xyz = solve(self.M_cat02, rgb)
        return dot(self.invM_, rgb_)


class CAM02(ColorSpace):
    def __init__(self, variant, c, Y_b, L_A, whitepoint=whitepoints_cie1931["D65"]):
        super().__init__(f"CAM02 ({variant})", ("J'", "a'", "b'"), 0)
        params = {
            "LCD": (0.77, 0.007, 0.0053),
            "SCD": (1.24, 0.007, 0.0363),
            "UCS": (1.00, 0.007, 0.0228),
        }
        self.K_L, self.c1, self.c2 = params[variant]
        self.ciecam02 = CIECAM02(c, Y_b, L_A, whitepoint)

    def from_xyz100(self, xyz):
        J, _, _, h, M, _, _ = self.ciecam02.from_xyz100(xyz)
        J_ = (1 + 100 * self.c1) * J / (1 + self.c1 * J)
        M_ = 1 / self.c2 * numpy.log(1 + self.c2 * M)
        h_ = h / 180 * numpy.pi
        return numpy.array([J_, M_ * numpy.cos(h_), M_ * numpy.sin(h_)])

    def to_xyz100(self, jab):
        J_, a, b = jab
        J = J_ / (1 - (J_ - 100) * self.c1)
        h = numpy.mod(numpy.arctan2(b, a), 2 * numpy.pi) / numpy.pi * 180
        M_ = numpy.hypot(a, b)
        M = (numpy.exp(M_ * self.c2) - 1) / self.c2
        return self.ciecam02.to_xyz100(numpy.array([J, M, h]), "JMh")
