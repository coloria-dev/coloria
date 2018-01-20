# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .illuminants import white_point, d65
from . import srgb_linear
from . import srgb1


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
    def __init__(self, c, Y_b, L_A, whitepoint=white_point(d65())):
        # step0: Calculate all values/parameters which are independent of input
        #        samples
        # Table 2.3 Surround parameters

        # Nc and F are modelled as a function of c, and can be linearly
        # interpolated.
        c_vals = [0.69, 0.59, 0.535]
        F_Nc_vals = [1.0, 0.9, 0.8]
        assert 0.535 <= c <= 0.69
        F = numpy.interp(c_vals, c, F_Nc_vals)
        self.c = c
        self.N_c = F

        XYZ_w = whitepoint
        self.M_cat02 = numpy.array([
            [+0.7328, +0.4296, -0.1624],
            [-0.7036, +1.6975, +0.0061],
            [+0.0030, +0.0136, +0.9834],
            ])
        rgb_w = numpy.dot(self.M_cat02, XYZ_w)

        D = F * (1 - 1/3.6 * numpy.exp((-L_A-42)/92))
        D = min(D, 1.0)
        D = max(D, 0.0)

        X_w, Y_w, Z_w = XYZ_w

        self.D_rgb = D*Y_w/rgb_w + 1 - D
        k = 1 / (5*L_A + 1)
        self.F_L = 0.2 * k**4 * 5*L_A + 0.1*(1-k**4)**2 * numpy.cbrt(5*L_A)

        n = Y_b / Y_w
        self.z = 1.48 + numpy.sqrt(n)
        self.N_bb = 0.725 * 1/n**0.2
        self.N_cb = self.N_bb

        rgb_wc = self.D_rgb * rgb_w

        self.M_hpe = numpy.array([
            [+0.38971, +0.68898, -0.07868],
            [-0.22981, +1.18340, +0.04641],
            [+0.00000, +0.00000, +1.00000],
            ])
        rgb_w_ = numpy.dot(
            self.M_hpe, numpy.linalg_solve(self.M_cat02, rgb_wc)
            )

        alpha = (self.F_L*rgb_w_/100)**0.42
        rgb_aw_ = 400 * (alpha / (alpha + 27.13)) + 0.1

        self.A_w = (numpy.dot([2, 1, 0.05], rgb_aw_) - 0.305) * self.N_bb

        self.hi = [20.14, 90.00, 164.25, 237.53, 380.14]
        self.ei = [0.8, 0.7, 1.0, 1.2, 0.8]
        self.Hi = [0.0, 100.0, 200.0, 300.0, 400.0]
        return

    def from_xyz(self, xyz):
        # TODO scale xyz w.r.t. test illuminant?

        # Step 1: Calculate (sharpened) cone responses (transfer
        #         colour-matching functions to sharper sensors)
        rgb = numpy.dot(self.M_cat02, xyz)

        # Step 2: Calculate the corresponding (sharpened) cone response
        #         (considering various luminance level and surround conditions
        #         included in D; hence, in DR, DG and DB)
        rgb_c = self.D_rgb * rgb

        # Step 3: Calculate the Hunt-Pointer-Estevez response
        rgb_ = numpy.dot(self.M_hpe, numpy.linalg.solve(self.M_cat02, rgb_c))

        # Step 4: Calculate the post-adaptation cone response (resulting in
        #         dynamic range compression)
        alpha = (self.F_L*abs(rgb_)/100)**0.42
        rgb_a_ = numpy.sign(rgb_[0]) * 400 * (alpha / (alpha+27.13)) + 0.1

        # Step 5: Calculate Redness–Greenness (a) , Yellowness–Blueness (b)
        #         components and hue angle (h)
        a = numpy.dot([1, -12/11, 1/11], rgb_a_)
        b = numpy.dot([1, 1, -2], rgb_a_) / 9
        h = numpy.arctan2(b/a)
        assert 0 < h < 2*numpy.pi

        # Step 6: Calculate eccentricity (et) and hue composition (H), using
        #         the unique hue data given in Table 2.4.
        # Red Yellow Green Blue Red
        #
        h_ = h + 2*numpy.pi if h < self.hi[0] else h
        e_t = 0.25 * (numpy.cos(h_+2) + 3.8)
        i0 = numpy.where(h_ > self.hi)[0][0] - 1
        beta = (h_ - self.hi[i0]) / self.ei[i0]
        H = self.Hi[i0] + 100 * beta \
            / (beta + (self.hi[i0+1] - h_)/self.ei[i0+1])

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

    def to_xyz(self, input_dict):
        # Step 1: Obtain J, C and h from H, Q, M, s
        #
        if 'J' in input_dict:
            J = input_dict['J']
            # Q perhaps needed for C
            Q = (4/self.c) * numpy.sqrt(J/100) * (self.A_w+4) * self.F_L**0.25
        else:
            # Step 1–1: Compute J from Q (if start from Q)
            Q = input_dict['Q']
            J = 6.25 * (self.c*Q / (self.A_w+4) / self.F_L**0.25)**2

        # Step 1–2: Calculate C from M or s
        if 'M' in input_dict:
            C = input_dict['M'] / self.F_L**0.25
        else:
            s = input_dict['s']
            C = (s/100)**2 * (Q/self.F_L**0.25)

        if 'h' in input_dict:
            h = input_dict['h']
        else:
            # Step 1–3: Calculate h from H (if start from H)
            H = input_dict['H']
            i0 = numpy.where(H > self.Hi)[0][0] - 1
            Hi = self.Hi[i0]
            hi, hi1 = self.hi[i0], self.hi[i0+1]
            ei, ei1 = self.ei[i0], self.ei[i0+1]
            h_ = ((H - Hi) * (ei1*hi - ei*hi1) - 100*hi*ei1) \
                / ((H - Hi) * (ei1 - ei) - 100*ei1)
            h = h_-2*numpy.pi if h_ > 2*numpy.pi else h_

        # Step 2: Calculate t, et , p1, p2 and p3
        t = (C / numpy.sqrt(J/100) / (1.64-0.29**self.n)**0.74)**(1/0.9)
        e_t = 0.25 * (numpy.cos(h+2) + 3.8)
        A = self.A_w * (J/100)**(1/self.c/self.z)

        p2 = A / self.N_bb + 0.305

        # Step 3: Calculate a and b
        if abs(t) < 1.0e-15:
            a = 0
            b = 0
        else:
            p1 = (50000/13 * self.N_c * self.N_cb) * e_t * (1/t)
            p3 = 21/20
            cosh = numpy.cos(h)
            sinh = numpy.sin(h)
            if abs(sinh) >= abs(cosh):
                p4 = p1 / sinh
                b = p2 * (2+p3) * 460 / \
                    (1403*p4 + (2+p3)*220*(cosh/sinh) - 27 + p3*6300)
                a = b * cosh/sinh
            else:
                p5 = p1 / cosh
                a = p2 * (2+p3) * 460 / \
                    (1403*p5 + (2+p3)*220 - (27 - p3*6300)*sinh/cosh)
                b = a * sinh/cosh

        # Step 4: Calculate RGB_a_
        rgb_a_ = numpy.dot(numpy.array([
            [460, 451, 288],
            [460, -891, -261],
            [460, -220, -6300]
            ]), numpy.array([p2, a, b])
            ) / 1403

        # Step 5: Calculate RGB_
        rgb_ = numpy.sign(rgb_a_ - 0.1) * 100/self.F_L * (
            (27.13 * abs(rgb_a_-0.1)) / (400 - abs(rgb_a_-0.1))
            )**(1/0.42)

        # Step 6: Calculate RC, GC and BC
        rgb_c = numpy.dot(self.M_cat02, numpy.linalg.solve(self.M_hpe, rgb_))

        # Step 7: Calculate R, G and B
        rgb = rgb_c / self.D_rgb

        # Step 8: Calculate X, Y and Z
        xyz = numpy.solve(self.M_cat02, rgb)
        # TODO scale xyz w.r.t. test illuminant?
        return xyz

    def srgb_gamut(self, filename='srgb-ciecam02.vtu', n=50):
        import meshio
        import meshzoo
        points, cells = meshzoo.cube(nx=n, ny=n, nz=n)
        pts = self.from_xyz(srgb_linear.to_xyz(points.T)).T
        rgb = srgb1.from_srgb_linear(points)
        meshio.write(
            filename,
            pts, {'tetra': cells},
            point_data={'srgb': rgb}
            )
        return
