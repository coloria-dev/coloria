"""
See

Ming Ronnier Luo,
CIE Chromatic Adaptation; Comparison of von Kries, CIELAB, CMCCAT97 and CAT02,
<https://doi.org/10.1007/978-3-642-27851-8_321-1>.
"""
import numpy as np
from numpy.typing import ArrayLike


def _von_kries_type(
    M: np.ndarray, whitepoint_source: ArrayLike, whitepoint_target: ArrayLike
):
    whitepoint_source = np.asarray(whitepoint_source)
    whitepoint_target = np.asarray(whitepoint_target)

    rgb_rw = M @ whitepoint_target
    rgb_w = M @ whitepoint_source
    d = rgb_rw / rgb_w

    A = np.linalg.solve(M, (M.T * d).T)
    Ainv = np.linalg.solve(M, (M.T / d).T)
    return A, Ainv


def von_kries(whitepoint_source: ArrayLike, whitepoint_target: ArrayLike):
    # Different matrices are presented for the "von-Kries" transform.
    #
    # <https://doi.org/10.1007/978-3-642-27851-8_321-1>:
    #   [0.0, 1.0, 0.0],
    #   [-0.46, 1.36, 0.1],
    #   [0.0, 0.0, 1.0],
    #
    # <http://www.brucelindbloom.com/index.html?Eqn_ChromAdapt.html>:
    #    0.4002400  0.7076000 -0.0808100
    #   -0.2263000  1.1653200  0.0457000
    #    0.0000000  0.0000000  0.9182200
    #
    # <https://doi.org/10.1117/12.410788>,
    # <http://dx.doi.org/10.1002/col.20573>
    # (Hunt-Pointer-Estevez matrix):
    #    0.3897 0.6890 -0.0787
    #   -0.2298 1.1834  0.0464
    #    0      0       1
    #
    judd = np.array(
        [
            [0.0, 1.0, 0.0],
            [-0.46, 1.36, 0.1],
            [0.0, 0.0, 1.0],
        ]
    )
    return _von_kries_type(judd, whitepoint_source, whitepoint_target)
