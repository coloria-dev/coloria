import npx
import numpy as np

from colorio.illuminants import whitepoints_cie1931


class CAM16Legacy:
    """
    Legacy CAM16 implementation for comparison purposes.
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
        F = np.interp(c, c_vals, F_Nc_vals)
        self.c = c
        self.N_c = F

        self.M16 = np.array(
            [
                [+0.401288, +0.650173, -0.051461],
                [-0.250268, +1.204414, +0.045854],
                [-0.002079, +0.048952, +0.953127],
            ]
        )
        # The standard actually recommends using this approximation as
        # inversion operation.
        approx_inv_M16 = np.array(
            [
                [+1.86206786, -1.01125463, +0.14918677],
                [+0.38752654, +0.62144744, -0.00897398],
                [-0.01584150, -0.03412294, +1.04996444],
            ]
        )
        self.invM16 = np.linalg.inv(self.M16) if exact_inversion else approx_inv_M16
        RGB_w = np.dot(self.M16, whitepoint)

        D = F * (1 - 1 / 3.6 * np.exp((-L_A - 42) / 92))
        D = np.clip(D, 0.0, 1.0)

        self.D_RGB = D * Y_w / RGB_w + 1 - D

        k = 1 / (5 * L_A + 1)
        self.F_L = k**4 * L_A + 0.1 * (1 - k**4) ** 2 * np.cbrt(5 * L_A)

        self.n = Y_b / Y_w
        self.z = 1.48 + np.sqrt(self.n)
        self.N_bb = 0.725 / self.n**0.2
        self.N_cb = self.N_bb

        RGB_wc = self.D_RGB * RGB_w
        alpha = (self.F_L * RGB_wc / 100) ** 0.42
        RGB_aw_ = 400 * alpha / (alpha + 27.13) + 0.1
        self.A_w = (np.dot([2, 1, 1 / 20], RGB_aw_) - 0.305) * self.N_bb

        self.h = np.array([20.14, 90.00, 164.25, 237.53, 380.14])
        self.e = np.array([0.8, 0.7, 1.0, 1.2, 0.8])
        self.H = np.array([0.0, 100.0, 200.0, 300.0, 400.0])

    def from_xyz100(self, xyz):
        # Step 1: Calculate 'cone' responses
        rgb = npx.dot(self.M16, xyz)
        # Step 2: Complete the color adaptation of the illuminant in
        #         the corresponding cone response space
        rgb_c = (rgb.T * self.D_RGB).T

        # Step 3: Calculate the post-adaptation cone response (resulting in
        #         dynamic range compression)
        alpha = (self.F_L * abs(rgb_c) / 100) ** 0.42
        rgb_a = np.sign(rgb_c) * 400 * alpha / (alpha + 27.13) + 0.1

        # Step 4
        a = npx.dot(np.array([1, -12 / 11, 1 / 11]), rgb_a)
        b = npx.dot(np.array([1 / 9, 1 / 9, -2 / 9]), rgb_a)
        # Make sure that h is in [0, 360]
        h = np.rad2deg(np.arctan2(b, a)) % 360

        # Step 5: Calculate eccentricity (e_t) and hue composition (H), using
        #         the unique hue data given in Table 2.4.
        h_ = (h - self.h[0]) % 360 + self.h[0]
        e_t = (np.cos(np.deg2rad(h_) + 2) + 3.8) / 4
        i = np.searchsorted(self.h, h_) - 1
        beta = (h_ - self.h[i]) * self.e[i + 1]
        H = self.H[i] + 100 * beta / (beta + self.e[i] * (self.h[i + 1] - h_))

        # Step 6
        A = (npx.dot(np.array([2, 1, 1 / 20]), rgb_a) - 0.305) * self.N_bb

        # Step 7: Calculate the correlate of lightness
        J = 100 * (A / self.A_w) ** (self.c * self.z)

        # Step 8: Calculate the correlate of brightness
        sqrt_J_100 = np.sqrt(J / 100)
        Q = (4 / self.c) * sqrt_J_100 * (self.A_w + 4) * self.F_L**0.25

        # Step 9: Calculate the correlates of chroma (C), colourfulness (M)
        #          and saturation (s)
        #
        t = (
            50000
            / 13
            * e_t
            * self.N_c
            * self.N_cb
            * np.sqrt(a**2 + b**2)
            / npx.dot(np.array([1, 1, 21 / 20]), rgb_a)
        )
        C = t**0.9 * (1.64 - 0.29**self.n) ** 0.73 * sqrt_J_100
        M = C * self.F_L**0.25
        s = 100 * np.sqrt(M / Q)

        return np.array([J, C, H, h, M, s, Q])

    def to_xyz100(self, data, description):
        """Input: J or Q; C, M or s; H or h"""
        if description[0] == "J":
            J = data[0]
            # Q perhaps needed for C
            Q = (4 / self.c) * np.sqrt(J / 100) * (self.A_w + 4) * self.F_L**0.25
        else:
            # Step 1–1: Compute J from Q (if start from Q)
            assert description[0] == "Q"
            Q = data[0]
            J = 6.25 * (self.c * Q / (self.A_w + 4) / self.F_L**0.25) ** 2

        # Step 1–2: Calculate C from M or s
        if description[1] == "C":
            C = data[1]
        elif description[1] == "M":
            M = data[1]
            C = M / self.F_L**0.25
        else:
            assert description[1] == "s"
            s = data[1]
            C = (s / 100) ** 2 * Q / self.F_L**0.25

        if description[2] == "h":
            h = data[2]
        else:
            assert description[2] == "H"
            # Step 1–3: Calculate h from H (if start from H)
            H = data[2]
            i = np.searchsorted(self.H, H) - 1
            Hi = self.H[i]
            hi, hi1 = self.h[i], self.h[i + 1]
            ei, ei1 = self.e[i], self.e[i + 1]
            h_ = ((H - Hi) * (ei1 * hi - ei * hi1) - 100 * hi * ei1) / (
                (H - Hi) * (ei1 - ei) - 100 * ei1
            )
            h = np.mod(h_, 360)

        h = np.deg2rad(h)

        # Step 2: Calculate t, et , p1, p2 and p3
        A = self.A_w * (J / 100) ** (1 / self.c / self.z)

        # Step 3: Calculate a and b
        t = (C / np.sqrt(J / 100) / (1.64 - 0.29**self.n) ** 0.73) ** (1 / 0.9)
        e_t = 0.25 * (np.cos(h + 2) + 3.8)

        one_over_t = 1 / t
        one_over_t = np.select([np.isnan(one_over_t), True], [np.inf, one_over_t])

        p1 = (50000.0 / 13) * self.N_c * self.N_cb * e_t * one_over_t
        p2 = A / self.N_bb + 0.305
        p3 = 21 / 20

        sin_h = np.sin(h)
        cos_h = np.cos(h)

        num = p2 * (2 + p3) * (460.0 / 1403)
        denom_part2 = (2 + p3) * (220.0 / 1403)
        denom_part3 = (-27.0 / 1403) + p3 * (6300.0 / 1403)

        a = np.empty_like(h)
        b = np.empty_like(h)

        small_cos = np.abs(sin_h) >= np.abs(cos_h)
        b[small_cos] = num[small_cos] / (
            p1[small_cos] / sin_h[small_cos]
            + (denom_part2 * cos_h[small_cos] / sin_h[small_cos])
            + denom_part3
        )
        a[small_cos] = b[small_cos] * cos_h[small_cos] / sin_h[small_cos]
        a[~small_cos] = num[~small_cos] / (
            p1[~small_cos] / cos_h[~small_cos]
            + denom_part2
            + (denom_part3 * sin_h[~small_cos] / cos_h[~small_cos])
        )
        b[~small_cos] = a[~small_cos] * sin_h[~small_cos] / cos_h[~small_cos]

        # Step 4: Calculate RGB_a_
        rgb_a_ = (
            npx.dot(
                np.array([[460, 451, 288], [460, -891, -261], [460, -220, -6300]]),
                np.array([p2, a, b]),
            )
            / 1403
        )

        # Step 5: Calculate RGB_
        rgb_c = (
            np.sign(rgb_a_ - 0.1)
            * 100
            / self.F_L
            * ((27.13 * abs(rgb_a_ - 0.1)) / (400 - abs(rgb_a_ - 0.1))) ** (1 / 0.42)
        )

        # Step 6: Calculate R, G and B
        rgb = (rgb_c.T / self.D_RGB).T

        # Step 7: Calculate X, Y and Z
        xyz = npx.dot(self.invM16, rgb)
        return xyz
