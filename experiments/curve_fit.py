# -*- coding: utf-8 -*-
#
from __future__ import print_function

import matplotlib.pyplot as plt
import numpy
import scipy.optimize

import colorio


def f(x, *p):
    return p[0] * numpy.exp(-(x-p[1])**2 / p[2]**2)


def f2(x, *p):
    return (
        + p[0] * numpy.exp(-(x-p[1])**2 / p[2]**2)
        + p[3] * numpy.exp(-(x-p[4])**2 / p[5]**2)
        )


def cv():
    lmbda, data = colorio.observers.cie_1931_2()

    # You have to scale the x-data; see
    # <https://github.com/scipy/scipy/issues/8369>.
    lmbda *= 1.0e5

    popt1, _ = scipy.optimize.curve_fit(
            f, lmbda, data[1], p0=[1, 0, 1],
            )
    print(popt1)
    popt2, _ = scipy.optimize.curve_fit(
            f2, lmbda, data[1], p0=numpy.concatenate([popt1, [1, 0, 1]]),
            )
    print(popt2)

    plt.plot(lmbda, data[1])
    plt.plot(lmbda, f(lmbda, *popt1), ':')
    plt.plot(lmbda, f2(lmbda, *popt2), '.')
    plt.show()
    return


if __name__ == '__main__':
    cv()
