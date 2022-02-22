import npx
import numpy as np
from numpy.typing import ArrayLike

from ..illuminants import whitepoints_cie1931
from ._color_space import ColorSpace
from ._helpers import register


class RLAB(ColorSpace):
    """
    Mark D. Fairchild,
    Refinement of the RLAB Color Space,
    1996,
    <https://doi.org/10.1002/(SICI)1520-6378(199610)21:5<338::AID-COL3>3.0.CO;2-Z>.

    See also
    Mark D. Fairchild,
    Color Appearance Models, Second Edition,
    <https://last.hit.bme.hu/download/firtha/video/Colorimetry/Fairchild_M._Color_appearance_models__2005.pdf>.

    and

    Mark D. Fairchild,
    RLAB: A Color Appearance Space for Color Reproduction,
    Device-Independent Color Imaging and Imaging Systems Integration,
    1993, Proceedings Volume 1909,
    <https://doi.org/10.1117/12.149061>

    for the original RLAB. (This implementation uses the "refined RLAB" scheme.)
    """

    name = "RLAB"
    labels = ("LR", "aR", "bR")
    k0 = 0

    def __init__(
        self,
        Y_n: float = 318.0,
        D: float = 0.0,
        whitepoint: ArrayLike = whitepoints_cie1931["D65"],
        sigma: float = 1.0 / 2.3,
    ):
        # One purpose of RLAB is to account for the adaptation in the human visual
        # system. That is, the visual system sees a red apple no matter if it is looked
        # at in bright daylight, at dawn, or in the light of a fire.  To achieve this,
        # RLAB first converts a tristimulus value XYZ into a tristimulus value under the
        # reference illuminant D65.

        # <https://en.wikipedia.org/wiki/LMS_color_space#Hunt,_RLAB>.
        # M_HPE, normalized to equal-energy illuminant:
        M = np.array(
            [[0.3897, 0.6890, -0.0787], [-0.2298, 1.1834, 0.0464], [0.0, 0.0, 1.0]]
        )

        self.whitepoint = np.asarray(whitepoint)

        lms_n = M @ self.whitepoint
        lms_e = 3.0 * lms_n / np.sum(lms_n)
        Yn3 = Y_n ** (1 / 3)
        p_lms = (1.0 + Yn3 + lms_e) / (1.0 + Yn3 + 1.0 / lms_e)
        a_lms = (p_lms + D * (1 - p_lms)) / lms_n

        # reference stimulus
        # The matrix R below can be computed explicitly with
        #
        #   wp = whitepoints_cie1931["D65"]
        #   Y_n = 318.0
        #   D = 0.0
        #   #
        #   lms_n = M @ wp
        #   lms_e = 3 * lms_n / np.sum(lms_n)
        #   Yn3 = Y_n ** (1 / 3)
        #   p_lms = (1 + Yn3 + lms_e) / (1 + Yn3 + 1 / lms_e)
        #   a_lms = (p_lms + D * (1 - p_lms)) / lms_n
        #   # D65-normalized Hunt-Pointer-Estevez matrix:
        #   M = np.array([
        #       [0.4002, 0.7076, -0.0808],
        #       [-0.2263, 1.1653, 0.0457],
        #       [0.0, 0.0, 0.9182]
        #   ])
        #   R = np.diag(1 / wp) @ np.linalg.inv(M) @ np.diag(1 / a_lms)
        #
        # In the remainder of model, however, one uses the equal-energy normalized
        # matrix. This seems to be a mistake in the model.
        # TODO M looks like it already is the inverse of another matrix given to
        # some digits. Find out which.
        S = np.array(
            [
                [1.03565229, -0.05781178, 0.02215949],
                [-0.0060366, 1.00482272, 0.00121388],
                [0.0, 0.0, 1.0],
            ]
        )
        # R = np.array(
        #     [[1.9569, -1.1882, 0.2313], [0.3612, 0.6388, 0.0], [0.0, 0.0, 1.0]]
        # )
        # R = S @ np.linalg.inv(M)

        self.sigma = sigma

        self.A = S @ np.linalg.solve(M, (a_lms * M.T).T)
        self.Ainv = np.linalg.solve(M, ((M @ np.linalg.inv(S)).T / a_lms).T)

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
    def from_xyz100(self, xyz: ArrayLike) -> np.ndarray:
        # First, the stimuli xyz are translated into reference stimuli xyz_ref to
        # account for the environment adaptation of the human visual system.
        xyz_ref = npx.dot(self.A, xyz)

        x_ref_s, y_ref_s, z_ref_s = xyz_ref**self.sigma

        # LR represents an achromatic response analogous to CIELAB L*. The red-green
        # chromatic response is given by a_R (analogous to CIELAB a*) and the
        # yellow–blue chromatic response is given by b_R (analogous to CIELAB b*).
        L_R = 100 * y_ref_s
        a_R = 430 * (x_ref_s - y_ref_s)
        b_R = 170 * (y_ref_s - z_ref_s)

        return np.array([L_R, a_R, b_R])

    def to_xyz100(self, lab: ArrayLike) -> np.ndarray:
        L_R, a_R, b_R = np.asarray(lab)

        y_ref_s = L_R / 100
        x_ref_s = a_R / 430 + y_ref_s
        z_ref_s = y_ref_s - b_R / 170

        xyz_ref = np.array([x_ref_s, y_ref_s, z_ref_s]) ** (1.0 / self.sigma)
        return npx.dot(self.Ainv, xyz_ref)


register("rlab", RLAB())
