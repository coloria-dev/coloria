import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf

import colorio


def plot_orig():
    obs = colorio.observers.cie_1931_2()
    plt.plot(obs.lmbda_nm, obs.data[0], label="x", color="r")
    plt.plot(obs.lmbda_nm, obs.data[1], label="y", color="g")
    plt.plot(obs.lmbda_nm, obs.data[2], label="z", color="b")


def plot_wss():
    """Wyman, Sloan, Shirley
    Simple Analytic Approximations to the CIE XYZ Color Matching Function
    """
    lmbda = np.arange(360, 831)

    obs = colorio.observers.cie_1931_2()

    def g(x, alpha, mu, sigma1, sigma2):
        sigma = np.full(x.shape, sigma1)
        sigma[x > mu] = sigma2
        return alpha * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    # better fit, but negative:
    # x_ = (
    #     +g(lmbda, 0.363, 440.8, 15.0, 50.0)
    #     + g(lmbda, 1.056, 599.4, 36.2, 31.2)
    #     + g(lmbda, -0.212, 493.9, 20.5, 25.6)
    # )
    # plt.plot(lmbda, x_ - obs.data[0], label="x")

    # original data:
    x2_ = (
        +g(lmbda, 0.362, 442.0, 1 / 0.0624, 1 / 0.0374)
        + g(lmbda, 1.056, 599.8, 1 / 0.0264, 1 / 0.0323)
        + g(lmbda, -0.065, 501.1, 1 / 0.0490, 1 / 0.0382)
    )
    plt.plot(lmbda, x2_ - obs.data[0], label="x2")

    print(np.min(x2_))
    assert np.all(x2_ > 0.0)

    y_ = g(lmbda, 0.821, 568.8, 46.9, 40.5) + g(lmbda, 0.286, 530.9, 16.3, 31.1)
    plt.plot(lmbda, y_ - obs.data[1], label="y")

    z_ = g(lmbda, 1.217, 437.0, 11.8, 36.0) + g(lmbda, 0.681, 459.0, 26.0, 13.8)
    plt.plot(lmbda, z_ - obs.data[2], label="z")

    plt.legend()

    # def fun(lmbda, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3):
    #     return g(lmbda, a1, b1, c1, d1) + g(lmbda, a2, b2, c2, d2)+ g(lmbda, a3, b3, c3, d3)

    # params0 = [0.362, 442.0, 16.0, 26.0, 1.056, 599.8, 37.0, 30.9, -0.065,
    #         501.1, 20.4, 26.18]
    # out, _ = curve_fit(fun, obs.lmbda_nm, obs.data[0], params0)
    # print(list(out))
    # sol[0]: 0.363, 440.8, 15.02, 49.99, 1.056, 599.4, 36.19, 31.18, -0.212, 493.9, 20.45, 25.58]
    # sol[1]: 0.821, 568.8, 46.89, 40.49, 0.286, 530.9, 16.31, 31.08]
    # sol[2]: 1.217, 437.0, 11.84, 35.99, 0.681, 459.0, 25.96, 13.79


def plot_new():
    def fun(lmbda, alpha, lmbda0, sigma, beta, alpha2, lmbda02, sigma2, beta2):
        return (
            alpha
            * np.exp(-((lmbda - lmbda0) ** 2) / 2 / sigma**2)
            * (1 + erf(beta * (lmbda - lmbda0) / np.sqrt(2) / sigma))
        ) + (
            alpha2
            * np.exp(-((lmbda - lmbda02) ** 2) / 2 / sigma2**2)
            * (1 + erf(beta2 * (lmbda - lmbda02) / np.sqrt(2) / sigma2))
        )

    params0 = [0.67, 450.0, 65.0, 2.0, 0.67, 520.0, 65.0, 2.0]
    # [0.7792439227085641, 526.7725443457584, 54.83668090879216, 1.3129513447747418]

    obs = colorio.observers.cie_1931_2()
    out, _ = curve_fit(fun, obs.lmbda_nm, obs.data[1], params0)
    print(list(out))

    lmbda = np.arange(360, 831)
    y_ = fun(lmbda, *out)

    plt.plot(lmbda, y_ - obs.data[1])


# plot_orig()
plot_wss()
# plot_new()
plt.show()
