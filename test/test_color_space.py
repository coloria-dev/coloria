import matplotlib.pyplot as plt
import numpy

import colorio


def test_visible_slice():
    xyy = colorio.XYY()
    xyy.save_visible_slice("xyy-visible-slice.png", 2, 0.4)
    plt.close()

    cielab = colorio.CIELAB()
    cielab.show_visible_slice(0, 50)
    plt.close()
    cielab.save_visible_slice("cielab-visible-slice.png", 0, 50)
    plt.close()

    L_A = 64 / numpy.pi / 5
    cam16 = colorio.CAM16UCS(0.69, 20, L_A)
    cam16.save_visible_slice("cam16ucs-visible-slice.png", 0, 50)
    plt.close()
    return


def test_macadams():
    # cielab = colorio.CIELAB()
    # cielab.show_macadams(0, 50)
    # cieluv = colorio.CIELUV()
    # cieluv.show_macadams(0, 50)
    jzazbz = colorio.JzAzBz()
    jzazbz.show_macadams(0, 0.5)
    # jzazbz.show_luo_rigg(0, 0.1)
    # xyy = colorio.XYY()
    # xyy.show_macadams(1.5, k0=2)
    #
    # L_A = 64 / numpy.pi / 5
    # cam02 = colorio.CAM02("UCS", 0.69, 20, L_A)
    # cam02.show_macadams(0, 50)
    #
    # L_A = 64 / numpy.pi / 5
    # cam16 = colorio.CAM16UCS(0.69, 20, L_A)
    # # cam16.show_macadams(0, 60)
    # cam16.show_luo_rigg(0, 60, ellipse_scaling=2.0)
    return


if __name__ == "__main__":
    test_visible_slice()