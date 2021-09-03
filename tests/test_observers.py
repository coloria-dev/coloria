import matplotlib.pyplot as plt
import pytest

import colorio
from colorio import SpectralData


@pytest.mark.parametrize(
    "observer", [colorio.observers.cie_1931_2(), colorio.observers.cie_1964_10()]
)
def test_observers(observer):
    # For plot colors, take the SRGB approximation of the color that the
    # observer would perceive if a light spectrum hits its eye that corresponds
    # to the sensitivity spectrum.
    colors = []
    for k in range(3):
        out = colorio.illuminants.spectrum_to_xyz100(
            SpectralData(observer.lmbda_nm, observer.data[k]), observer=observer
        )
        out *= 100 / out[1]
        srgb = colorio.cs.SrgbLinear()
        rgb_vals = srgb.from_xyz100(out)
        rgb_vals[rgb_vals < 0] = 0
        # project down to proper rgb
        rgb_vals /= max(rgb_vals)
        colors.append(srgb.to_rgb1(rgb_vals))

    plt.plot(observer.lmbda_nm, observer.data[0], color=colors[0])
    plt.plot(observer.lmbda_nm, observer.data[1], color=colors[1])
    plt.plot(observer.lmbda_nm, observer.data[2], color=colors[2])

    plt.xlabel("wavelength [nm]")
    plt.grid()
    plt.xlim(observer.lmbda_nm[0], observer.lmbda_nm[-1])
    plt.ylim(ymin=0)

    plt.show()


if __name__ == "__main__":
    test_observers(
        # colorio.observers.cie_1931_2()
        colorio.observers.cie_1964_10()
    )
