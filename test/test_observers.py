# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt

import colorio

import numpy


def test_observers():
    lmbda, data = colorio.observers.cie_1931_2()

    # For plot colors, take the SRGB approximation of the color that the
    # observer would perceive if a light spectrum hits its eye that corresponds
    # to the sensitivity spectrum.
    colors = []
    for k in range(3):
        out = colorio.illuminants.spectrum_to_xyz100(
                (lmbda, data[k]), observer=colorio.observers.cie_1931_2()
                )
        out *= 100 / out[1]
        srgb = colorio.SrgbLinear()
        rgb_vals = srgb.from_xyz100(out)
        rgb_vals[rgb_vals < 0] = 0
        # project down to proper rgb
        rgb_vals /= max(rgb_vals)
        colors.append(srgb.to_srgb1(rgb_vals))

    plt.plot(lmbda, data[0], color=colors[0])
    plt.plot(lmbda, data[1], color=colors[1])
    plt.plot(lmbda, data[2], color=colors[2])

    plt.xlabel('wavelength (nm)')
    plt.grid()
    plt.xlim(lmbda[0], lmbda[-1])
    plt.ylim(ymin=0)

    plt.show()
    return


if __name__ == '__main__':
    test_observers()
