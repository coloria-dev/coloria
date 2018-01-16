# -*- coding: utf-8 -*-
#
import colorio

import matplotlib.pyplot as plt


def test_observers():
    lmbda, data = colorio.observers.cie_standard_observer_2()

    plt.plot(lmbda, data[0])
    plt.plot(lmbda, data[1])
    plt.plot(lmbda, data[2])

    plt.xlabel('wavelength (nm)')
    plt.grid()
    plt.xlim(lmbda[0], lmbda[-1])
    plt.ylim(ymin=0)

    plt.show()
    return


if __name__ == '__main__':
    test_observers()
