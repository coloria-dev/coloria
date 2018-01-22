# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .illuminants import white_point, d65


def find_first(a, alpha):
    '''Given an array a and a value alpha, this method returns the first index
    i where a[i] > alpha. Vectorized in alpha.
    '''
    # https://stackoverflow.com/a/48367770/353337
    return numpy.argmax(numpy.add.outer(alpha, -a) < 0, axis=-1)


class CIECAM02(object):
    '''
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
    '''
    # pylint: disable=too-many-instance-attributes
    def __init__(self, c, Y_b, L_A, whitepoint=white_point(d65())):
        # step0: Calculate all values/parameters which are independent of input
        #        samples
        Y_w = whitepoint[1]

        # Nc and F are modelled as a function of c, and can be linearly
        # interpolated.
        c_vals = [0.535, 0.59, 0.69]
        F_Nc_vals = [0.8, 0.9, 1.0]
        assert 0.535 <= c <= 0.69
        F = numpy.interp(c, c_vals, F_Nc_vals)
        self.c = c
        self.N_c = F

        self.M_cat02 = numpy.array([
            [+0.7328, +0.4296, -0.1624],
            [-0.7036, +1.6975, +0.0061],
            [+0.0030, +0.0136, +0.9834],
            ])
        RGB_w = numpy.dot(self.M_cat02, whitepoint)

        D = F * (1 - 1/3.6 * numpy.exp((-L_A-42)/92))
        D = min(D, 1.0)
        D = max(D, 0.0)

        self.D_RGB = D*Y_w/RGB_w + 1 - D

        k = 1 / (5*L_A + 1)
        self.F_L = k**4 * L_A + 0.1*(1-k**4)**2 * numpy.cbrt(5*L_A)

        self.n = Y_b / Y_w
        self.z = 1.48 + numpy.sqrt(self.n)
        self.N_bb = 0.725 / self.n**0.2
        self.N_cb = self.N_bb

        RGB_wc = self.D_RGB * RGB_w

        self.M_hpe = numpy.array([
            [+0.38971, +0.68898, -0.07868],
            [-0.22981, +1.18340, +0.04641],
            [+0.00000, +0.00000, +1.00000],
            ])
        RGB_w_ = numpy.dot(
            self.M_hpe, numpy.linalg.solve(self.M_cat02, RGB_wc)
            )

        alpha = (self.F_L*RGB_w_/100)**0.42
        RGB_aw_ = 400 * alpha / (alpha + 27.13) + 0.1

        self.A_w = (numpy.dot([2, 1, 1/20], RGB_aw_) - 0.305) * self.N_bb

        self.h = \
            numpy.array([20.14, 90.00, 164.25, 237.53, 380.14]) * numpy.pi/180
        self.e = numpy.array([0.8, 0.7, 1.0, 1.2, 0.8])
        self.H = numpy.array([0.0, 100.0, 200.0, 300.0, 400.0])
        return

    def from_xyz(self, xyz):
        # TODO scale xyz w.r.t. test illuminant?

        # Step 1: Calculate (sharpened) cone responses (transfer
        #         colour-matching functions to sharper sensors)
        rgb = numpy.dot(self.M_cat02, xyz)

        # Step 2: Calculate the corresponding (sharpened) cone response
        #         (considering various luminance level and surround conditions
        #         included in D; hence, in DR, DG and DB)
        rgb_c = (rgb.T * self.D_RGB).T

        # Step 3: Calculate the Hunt-Pointer-Estevez response
        rgb_ = numpy.dot(self.M_hpe, numpy.linalg.solve(self.M_cat02, rgb_c))

        # Step 4: Calculate the post-adaptation cone response (resulting in
        #         dynamic range compression)
        alpha = (self.F_L*abs(rgb_)/100)**0.42
        rgb_a_ = numpy.sign(rgb_) * 400 * alpha / (alpha+27.13) + 0.1

        # Step 5: Calculate Redness–Greenness (a) , Yellowness–Blueness (b)
        #         components and hue angle (h)
        a = numpy.dot([1, -12/11, 1/11], rgb_a_)
        b = numpy.dot([1, 1, -2], rgb_a_) / 9
        # Make sure that h is in [0, 2*pi]
        h = numpy.mod(numpy.arctan2(b, a), 2*numpy.pi)
        assert numpy.all(h >= 0) and numpy.all(h <= 2*numpy.pi)

        # Step 6: Calculate eccentricity (e_t) and hue composition (H), using
        #         the unique hue data given in Table 2.4.
        h_ = numpy.mod(h - self.h[0], 2*numpy.pi) + self.h[0]
        assert numpy.all(self.h[0] <= h_) and numpy.all(h_ < self.h[-1])
        e_t = 1/4 * (numpy.cos(h_+2) + 3.8)
        i = find_first(self.h, h_) - 1
        assert numpy.all(self.h[i] <= h_) and numpy.all(h_ <= self.h[i+1])
        beta = (h_ - self.h[i]) / self.e[i]
        H = self.H[i] + 100 * beta / (beta + (self.h[i+1] - h_)/self.e[i+1])

        # Step 7: Calculate achromatic response A
        A = (numpy.dot([2, 1, 1/20], rgb_a_) - 0.305) * self.N_bb

        # Step 8: Calculate the correlate of lightness
        J = 100 * (A/self.A_w)**(self.c*self.z)

        # Step 9: Calculate the correlate of brightness
        Q = (4/self.c) * numpy.sqrt(J/100) * (self.A_w + 4) * self.F_L**0.25

        # Step 10: Calculate the correlates of chroma (C), colourfulness (M)
        #          and saturation (s)
        t = (50000/13 * self.N_c * self.N_cb) * e_t * numpy.sqrt(a**2 + b**2) \
            / numpy.dot([1, 1, 21/20], rgb_a_)
        C = t**0.9 * numpy.sqrt(J/100) * (1.64 - 0.29**self.n)**0.73
        M = C * self.F_L**0.25
        s = 100 * numpy.sqrt(M/Q)
        return numpy.array([J, C, H, h, M, s, Q])

    def to_xyz(self, data, description):
        '''Input: J or Q; C, M or s; H or h
        '''
        # Step 1: Obtain J, C and h from H, Q, M, s
        #
        if description[0] == 'J':
            J = data[0]
            # Q perhaps needed for C
            Q = (4/self.c) * numpy.sqrt(J/100) * (self.A_w+4) * self.F_L**0.25
        else:
            # Step 1–1: Compute J from Q (if start from Q)
            assert description[0] == 'Q'
            Q = data[0]
            J = 6.25 * (self.c*Q / (self.A_w+4) / self.F_L**0.25)**2

        # Step 1–2: Calculate C from M or s
        if description[1] == 'C':
            C = data[1]
        elif description[1] == 'M':
            M = data[1]
            C = M / self.F_L**0.25
        else:
            assert description[1] == 's'
            s = data[1]
            C = (s/100)**2 * Q / self.F_L**0.25

        if description[2] == 'h':
            h = data[2]
        else:
            assert description[2] == 'H'
            # Step 1–3: Calculate h from H (if start from H)
            H = data[2]
            i = find_first(self.H, H) - 1
            assert numpy.all(self.H[i] <= H) and numpy.all(H < self.H[i+1])
            Hi = self.H[i]
            hi, hi1 = self.h[i], self.h[i+1]
            ei, ei1 = self.e[i], self.e[i+1]
            h_ = ((H - Hi) * (ei1*hi - ei*hi1) - 100*hi*ei1) \
                / ((H - Hi) * (ei1 - ei) - 100*ei1)
            h = numpy.mod(h_, 2*numpy.pi)

        # Step 2: Calculate t, et , p1, p2 and p3
        t = (C / numpy.sqrt(J/100) / (1.64-0.29**self.n)**0.73)**(1/0.9)
        e_t = 0.25 * (numpy.cos(h+2) + 3.8)
        A = self.A_w * (J/100)**(1/self.c/self.z)

        p2 = A / self.N_bb + 0.305

        # Step 3: Calculate a and b
        # ENH In the specification, the case t=0 is treated separately as
        # otherwise, division by 0 occurs in
        #     p1 = 50000/13 * self.N_c * self.N_cb * e_t / t.
        # It turns out that things simplify a great deal when simply
        # calculating with one_over_p1. (This is legal since all other
        # quantities in the above term are nonzero.)
        one_over_p1 = t / e_t * 13/50000 / self.N_c / self.N_cb
        p3 = 21/20
        cosh = numpy.cos(h)
        sinh = numpy.sin(h)
        one_over_p4 = one_over_p1 * sinh
        one_over_p5 = one_over_p1 * cosh
        a, b = numpy.array([one_over_p5, one_over_p4]) * (
            p2 * (2+p3) * 460 /
            (1403 + one_over_p5*(2+p3)*220 - one_over_p4*(27 - p3*6300))
            )

        # Step 4: Calculate RGB_a_
        rgb_a_ = numpy.dot(numpy.array([
            [460, 451, 288],
            [460, -891, -261],
            [460, -220, -6300]
            ]), numpy.array([p2, a, b])) / 1403

        # Step 5: Calculate RGB_
        rgb_ = numpy.sign(rgb_a_ - 0.1) * 100/self.F_L * (
            (27.13 * abs(rgb_a_-0.1)) / (400 - abs(rgb_a_-0.1))
            )**(1/0.42)

        # Step 6: Calculate RC, GC and BC
        rgb_c = numpy.dot(self.M_cat02, numpy.linalg.solve(self.M_hpe, rgb_))

        # Step 7: Calculate R, G and B
        rgb = (rgb_c.T / self.D_RGB).T

        # Step 8: Calculate X, Y and Z
        xyz = numpy.linalg.solve(self.M_cat02, rgb)
        # TODO scale xyz w.r.t. test illuminant?
        return xyz


class CAM02(object):
    # pylint: disable=too-many-arguments
    def __init__(self, variant, c, Y_b, L_A, whitepoint=white_point(d65())):
        params = {
            'LCD': (0.77, 0.007, 0.0053),
            'SCD': (1.24, 0.007, 0.0363),
            'UCS': (1.00, 0.007, 0.0228),
            }
        self.K_L, self.c1, self.c2 = params[variant]
        self.ciecam02 = CIECAM02(c, Y_b, L_A, whitepoint)
        return

    def from_xyz(self, xyz):
        J, _, _, h, M, _, _ = self.ciecam02.from_xyz(xyz)
        J_ = (1+100*self.c1)*J / (1 + self.c1*J)
        M_ = 1/self.c2 * numpy.log(1 + self.c2*M)
        return numpy.array([J_, M_*numpy.cos(h), M_*numpy.sin(h)])

    def to_xyz(self, jab):
        J_, a, b = jab
        J = J_ / (1 - (J_-100)*self.c1)
        h = numpy.arctan2(b, a)
        M_ = numpy.sqrt(a**2 + b**2)
        M = (numpy.exp(M_ * self.c2) - 1) / self.c2
        return self.ciecam02.to_xyz(numpy.array([J, M, h]), 'JMh')
