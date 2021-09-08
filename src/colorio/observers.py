import json
import pathlib

import numpy as np

from ._helpers import SpectralData

this_dir = pathlib.Path(__file__).resolve().parent


def cie_1931_2(stepsize: int = 1):
    return _from_file(this_dir / "data/observers/cie-1931-2.json", stepsize)


def cie_1964_10(stepsize: int = 1):
    return _from_file(this_dir / "data/observers/cie-1964-10.json", stepsize)


def _from_file(filename: pathlib.Path, stepsize: int):
    with open(filename) as f:
        data = json.load(f)

    lmbda_start, lmbda_end, lmbda_step = data["lambda_nm"]
    assert lmbda_step == 1
    lmbda = np.arange(lmbda_start, lmbda_end + 1, stepsize)
    vals = np.array(data["xyz"])[:, ::stepsize]

    return SpectralData(lmbda, vals, data["name"])


def wws_cie_1931_2(lmbda):
    """
    Wyman, Sloan, Shirley,
    Simple Analytic Approximations to the CIE XYZ Color Matching Function,
    Journal of Computer Graphics Techniques,
    Vol. 2, No. 2, 2013,
    <http://cwyman.org/papers/jcgt13_xyzApprox.pdf>.
    """

    def g(x, alpha, mu, sigma1, sigma2):
        sigma = np.full(x.shape, sigma1)
        sigma[x > mu] = sigma2
        return alpha * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    # better fit, but has negative values in the range:
    # x_ = (
    #     +g(lmbda, 0.363, 440.8, 15.0, 50.0)
    #     + g(lmbda, 1.056, 599.4, 36.2, 31.2)
    #     + g(lmbda, -0.212, 493.9, 20.5, 25.6)
    # )
    # plt.plot(lmbda, x_ - obs.data[0], label="x")

    # original data:
    x_ = (
        g(lmbda, 0.362, 442.0, 16.0, 26.7)
        + g(lmbda, 1.056, 599.8, 37.9, 31.0)
        + g(lmbda, -0.065, 501.1, 20.4, 26.2)
    )

    y_ = g(lmbda, 0.821, 568.8, 46.9, 40.5) + g(lmbda, 0.286, 530.9, 16.3, 31.1)

    z_ = g(lmbda, 1.217, 437.0, 11.8, 36.0) + g(lmbda, 0.681, 459.0, 26.0, 13.8)

    return np.array([x_, y_, z_])


def wws_cie_1964_10(lmbda):
    x_ = +0.398 * np.exp(-1250 * np.log((lmbda + 570.1) / 1014) ** 2) + 1.132 * np.exp(
        -234 * np.log((1338 - lmbda) / 734.5) ** 2
    )
    y_ = 1.011 * np.exp(-0.5 * ((lmbda - 556.1) / 46.14) ** 2)
    z_ = 2.060 * np.exp(-32.0 * np.log((lmbda - 265.8) / 180.4) ** 2)
    return np.array([x_, y_, z_])
