import numpy


# parameters:
# graphic arts:
#   k_L = 1
#   K_1 = 0.045
#   K_2 = 0.015
#
# textiles:
#   k_L = 2
#   K_1 = 0.048
#   K_2 = 0.014
def cie94(lab1, lab2, k_L=1.0, K_1=0.045, K_2=0.015):
    l1, a1, b1 = numpy.asarray(lab1)
    l2, a2, b2 = numpy.asarray(lab2)

    dL = l1 - l2
    C1 = numpy.sqrt(a1 ** 2 + b1 ** 2)
    C2 = numpy.sqrt(a2 ** 2 + b2 ** 2)
    dC = C1 - C2
    da = a1 - a2
    db = b1 - b2
    # dH2 is mathematically >=0, but round-off can lead to small negatives
    dH2 = da ** 2 + db ** 2 - dC ** 2
    dH2 = numpy.maximum(dH2, numpy.zeros_like(dH2))
    dH = numpy.sqrt(dH2)

    S_L = 1.0
    S_C = 1 + K_1 * C1
    S_H = 1 + K_2 * C1

    k_C = 1.0
    k_H = 1.0
    dE = numpy.sqrt(
        (dL / k_L / S_L) ** 2 + (dC / k_C / S_C) ** 2 + (dH / k_H / S_H) ** 2
    )
    return dE
