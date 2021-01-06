import numpy

from .._linalg import dot, solve
from ..illuminants import whitepoints_cie1931
from ._color_space import ColorSpace


class RLAB(ColorSpace):
    """
    Mark D. Fairchild,
    Refinement of the RLAB Color Space,
    1996,
    <https://doi.org/10.1002/(SICI)1520-6378(199610)21:5<338::AID-COL3>3.0.CO;2-Z>.

    See also
    Mark D. Fairchild,
    Color Appearance Models, Second Edition,
    <http://last.hit.bme.hu/download/firtha/video/Colorimetry/Fairchild_M._Color_appearance_models__2005.pdf>.

    and

    Mark D. Fairchild,
    RLAB: A Color Appearance Space for Color Reproduction,
    Device-Independent Color Imaging and Imaging Systems Integration,
    1993, Proceedings Volume 1909,
    <https://doi.org/10.1117/12.149061>

    for the original RLAB. (This implementation uses the "refined RLAB" scheme.)
    """

    def __init__(
        self, Y_n=318.0, D=0.0, whitepoint=whitepoints_cie1931["D65"], sigma=1.0 / 2.3
    ):
        super().__init__("RLAB", ("LR", "aR", "bR"), 0)
        # One purpose of RLAB is to account for the adaptation in the human visual
        # system. That is, the visual system sees a red apple no matter if it is looked
        # at in bright daylight, at dawn, or in the light of a fire.  To achieve this,
        # RLAB first converts a tristimulus value XYZ into a tristimulus value under the
        # reference illuminant D65.

        # M_HPE, normalized to D65:
        # M_HPE, normalized to equal-energy illuminant,
        # <https://en.wikipedia.org/wiki/LMS_color_space#Hunt,_RLAB>.
        self.M = numpy.array(
            [[0.3897, 0.6890, -0.0787], [-0.2298, 1.1834, 0.0464], [0.0, 0.0, 1.0]]
        )

        lms_n = self.M @ whitepoint
        lms_e = 3.0 * lms_n / numpy.sum(lms_n)
        Yn3 = Y_n ** (1 / 3)
        p_lms = (1.0 + Yn3 + lms_e) / (1.0 + Yn3 + 1.0 / lms_e)
        self.a_lms = (p_lms + D * (1 - p_lms)) / lms_n

        # reference stimul
        # The matrix R below can be computed explicitly with
        #
        #   wp = whitepoints_cie1931["D65"]
        #   Y_n = 318.0
        #   D = 0.0
        #   #
        #   lms_n = self.M @ wp
        #   lms_e = 3 * lms_n / numpy.sum(lms_n)
        #   Yn3 = Y_n ** (1 / 3)
        #   p_lms = (1 + Yn3 + lms_e) / (1 + Yn3 + 1 / lms_e)
        #   a_lms = (p_lms + D * (1 - p_lms)) / lms_n
        #   self.R = numpy.diag(1 / wp) @ numpy.linalg.inv(M) @ numpy.diag(1 / a_lms)
        #
        # with M being the D65-normalized Hunt-Pointer-Estevez matrix.
        #
        # self.M = numpy.array([
        #     [0.4002, 0.7076, -0.0808],
        #     [-0.2263, 1.1653, 0.0457],
        #     [0.0, 0.0, 0.9182]
        # ]).
        #
        # In the remainder of model, however, one uses the equal-energy normalized
        # matrixi. This seems to be a mistake in the model.
        # TODO self.M looks like it already is the inverse of another matrix given to
        # some digits. Find out which.
        self.R = numpy.array(
            [[1.9569, -1.1882, 0.2313], [0.3612, 0.6388, 0.0], [0.0, 0.0, 1.0]]
        )

        self.sigma = sigma

    # Y_n is the absolute luminance of the adapting stimulus (typically a stimulus that
    # appears white in the image) in cd/m2 (=nit). Most laptop displays have a
    # brightness of 300 to 400 nits (2019).
    #
    # The D factor allows various proportions of cognitive discounting-the-illuminant. D
    # should be set equal to 1.0 for hard-copy images, 0.0 for soft-copy displays, and
    # an intermediate value for situations such as projected transparencies in
    # completely darkened rooms.
    #
    # The exponents in Equations 13.15–13.17 vary, depending on the relative luminance
    # of the surround. For an average surround sigma = 1/2.3, for a dim surround sigma =
    # 1/2.9, and for a dark surround sigma = 1/3.5.
    def from_xyz100(self, xyz):
        # First, the stimuli xyz are translated into reference stimuli xyz_ref to
        # account for the environment adaptation of the human visual system.
        lms_dash = (self.a_lms * dot(self.M, xyz).T).T
        xyz_ref = dot(self.R, lms_dash)

        x_ref_s, y_ref_s, z_ref_s = xyz_ref ** self.sigma

        # LR represents an achromatic response analogous to CIELAB L*. The redgreen
        # chromatic response is given by a R (analogous to CIELAB a*) and the
        # yellow–blue chromatic response is given by bR (analogous to CIELAB b*).
        L_R = 100 * y_ref_s
        a_R = 430 * (x_ref_s - y_ref_s)
        b_R = 170 * (y_ref_s - z_ref_s)

        return numpy.array([L_R, a_R, b_R])

    def to_xyz100(self, lab):
        L_R, a_R, b_R = lab

        y_ref_s = L_R / 100
        x_ref_s = a_R / 430 + y_ref_s
        z_ref_s = y_ref_s - b_R / 170

        xyz_ref = numpy.array([x_ref_s, y_ref_s, z_ref_s]) ** (1.0 / self.sigma)
        return solve(self.M, (solve(self.R, xyz_ref).T / self.a_lms).T)
