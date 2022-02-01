"""
Curve-fit observer curves with sums of exp(-x**2).
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import scipy.special

import colorio


def f(x, *p):
    return p[0] / np.sqrt(2 * np.pi) * np.exp(-((x - p[1]) ** 2) / (2 * p[2] ** 2))


def f3(x, *p):
    return (
        +p[0] / np.sqrt(2 * np.pi) * np.exp(-((x - p[1]) ** 2) / (2 * p[2] ** 2))
        + p[3] / np.sqrt(2 * np.pi) * np.exp(-((x - p[4]) ** 2) / (2 * p[5] ** 2))
        + p[6] / np.sqrt(2 * np.pi) * np.exp(-((x - p[7]) ** 2) / (2 * p[8] ** 2))
    )


# With a skew <https://en.wikipedia.org/wiki/Skew_normal_distribution>
def g(x, *p):
    return (
        p[0]
        * 1
        / np.sqrt(2 * np.pi)
        * np.exp(-((x - p[1]) ** 2) / (2 * p[2] ** 2))
        * (1 + scipy.special.erf(p[3] * (x - p[1]) / 2 / p[2]))
    )


def pade(x, *p):
    return (p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3) / (
        p[4] + p[5] * x + p[6] * x**2 + p[7] * x**3
    )


def cv():
    lmbda, data = colorio.observers.cie_1931_2()
    data = data[1]

    # You have to scale the x-data; see
    # <https://github.com/scipy/scipy/issues/8369>.
    lmbda *= 1.0e5

    popt1, _ = scipy.optimize.curve_fit(f, lmbda, data, p0=[1, 0, 1])
    print(popt1)

    popt2, _ = scipy.optimize.curve_fit(g, lmbda, data, p0=[1, 0, 1, 0], maxfev=10000)
    print(popt2)

    popt3, _ = scipy.optimize.curve_fit(
        pade, lmbda, data, p0=[1, 0, 0, 0, 1, 0, 0, 0], maxfev=10000
    )
    print(popt3)

    popt4, _ = scipy.optimize.curve_fit(
        f3, lmbda, data, p0=[1, 0, 1, 1, 0, 1, 1, 0, 1], maxfev=10000
    )
    print(popt4)

    plt.plot(lmbda, data, ".")
    plt.plot(lmbda, f(lmbda, *popt1), ":")
    plt.plot(lmbda, g(lmbda, *popt2), ":")
    plt.plot(lmbda, pade(lmbda, *popt3), ":")
    plt.plot(lmbda, f3(lmbda, *popt4), ":")
    plt.show()
    return


if __name__ == "__main__":
    cv()
