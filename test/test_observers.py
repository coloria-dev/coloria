# -*- coding: utf-8 -*-
#
import colorio

import matplotlib.pyplot as plt


def test_observers():
    lmbda, data = colorio.observers.cie_standard_observer_2()

    # Choose the colors such that they roughly correspond to the color of the
    # maximum response,
    # data[0]: 595 nm
    # data[1]: 557 nm
    # data[2]: 445 nm
    #
    # TODO adjust this better
    plt.plot(lmbda, data[0], color='#FFC400')
    plt.plot(lmbda, data[1], color='#ABFF00')
    plt.plot(lmbda, data[2], color='#0019FF')

    plt.xlabel('wavelength (nm)')
    plt.grid()
    plt.xlim(lmbda[0], lmbda[-1])
    plt.ylim(ymin=0)

    plt.show()
    return


if __name__ == '__main__':
    test_observers()
