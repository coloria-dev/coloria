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


def test_macadam():
    xyy = colorio.XYY()
    xyy.show_macadam(2, 0.4)
    plt.close()
    xyy.save_macadam("macadam-xyy.png", 2, 0.4)
    plt.close()

    cieluv = colorio.CIELUV()
    cieluv.save_macadam("macadam-cieluv.png", 0, 50)
    plt.close()

    jzazbz = colorio.JzAzBz()
    jzazbz.save_macadam("macadam-jzazbz.png", 0, 0.5)
    plt.close()

    # jzazbz.show_luo_rigg(0, 0.1)
    #
    # L_A = 64 / numpy.pi / 5
    # cam02 = colorio.CAM02("UCS", 0.69, 20, L_A)
    # cam02.show_macadam(0, 50)
    #
    # L_A = 64 / numpy.pi / 5
    # cam16 = colorio.CAM16UCS(0.69, 20, L_A)
    # # cam16.show_macadam(0, 60)
    # cam16.show_luo_rigg(0, 60, ellipse_scaling=2.0)
    return


def test_luo_rigg():
    xyy = colorio.XYY()
    xyy.show_luo_rigg(2, 0.4)
    plt.close()
    xyy.save_luo_rigg("luo-rigg-xyy.png", 2, 0.4)
    plt.close()

    cieluv = colorio.CIELUV()
    cieluv.save_luo_rigg("luo-rigg-cieluv.png", 0, 50, plot_srgb_gamut=False)
    plt.close()

    jzazbz = colorio.JzAzBz()
    jzazbz.save_luo_rigg("luo-rigg-jzazbz.png", 0, 0.5)
    plt.close()
    return


if __name__ == "__main__":
    # test_visible_slice()
    # test_macadam()
    test_luo_rigg()
