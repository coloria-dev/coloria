import numpy

from ..cs import CIELAB, CIELCH


# Acceptability: l=2, c=1 Perceptability: l=1, c=1
def cmc(lab1, lab2, l=2.0, c=1.0):
    cielab = CIELAB()
    cielch = CIELCH()
    lch1 = cielch.from_xyz100(cielab.to_xyz100(lab1))
    lch2 = cielch.from_xyz100(cielab.to_xyz100(lab2))

    L1, C1, h1 = lch1
    L2, C2, h2 = lch2

    F = numpy.sqrt(C1 ** 4 / (C1 ** 4 + 1900))
    if h1 >= 164 and h1 <= 345:
        T = 0.56 + numpy.abs(0.2 * numpy.cos(h1 + 168))
    else:
        T = 0.36 + numpy.abs(0.4 * numpy.cos(h1 + 35))

    if L1 < 16:
        S_L = 0.511
    else:
        S_L = (0.04095 * L1) / (1 + 0.01765 * L1)

    S_C = 0.0638 * C1 / (1 + 0.0131 * C1) + 0.638
    S_H = S_C * (F * T + 1 - F)

    dE = numpy.sqrt(
        ((L2 - L1) / l / S_L) ** 2 + ((C2 - C2) / c / S_C) ** 2 + ((h2 - h1) / S_H) ** 2
    )
    return dE
