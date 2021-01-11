import numpy


def ciede2000(lab1, lab2, k_L=1.0, k_C=1.0, k_H=1.0):
    L1, a1, b1 = lab1
    L2, a2, b2 = lab2

    dLp = L1 - L2

    Lpbar = (L1 + L2) / 2

    C1 = numpy.sqrt(a1 ** 2 + b1 ** 2)
    C2 = numpy.sqrt(a2 ** 2 + b2 ** 2)

    Cbar = (C1 + C2) / 2

    a1p = a1 + a1 / 2 * (1 - numpy.sqrt(Cbar ** 7 / (Cbar ** 7 + 25 ** 7)))
    a2p = a2 + a2 / 2 * (1 - numpy.sqrt(Cbar ** 7 / (Cbar ** 7 + 25 ** 7)))

    Cp1 = numpy.sqrt(a1p ** 2 + b1 ** 2)
    Cp2 = numpy.sqrt(a2p ** 2 + b2 ** 2)
    Cpbar = (Cp1 + Cp2) / 2
    dCp = C2 - C1

    hp1 = numpy.mod(numpy.arctan2(b1, a1p) / numpy.pi * 180, 360)
    hp2 = numpy.mod(numpy.arctan2(b2, a2p) / numpy.pi * 180, 360)

    dhp = numpy.abs(hp1 - hp2)

    if dhp <= 180:
        dhp = hp2 - hp1
    else:
        if hp2 <= hp1:
            dhp = hp2 - hp1 + 360
        else:
            dhp = hp2 - hp1 - 360

    dHp = 2 * numpy.sqrt(Cp1 * Cp2) * numpy.sin(dhp / 2)
    if numpy.abs(hp1 - hp2) <= 180:
        Hpbar = (hp1 + hp2) / 2
    else:
        if hp1 + hp2 < 360:
            Hpbar = (hp1 + hp2 + 360) / 2
        else:
            Hpbar = (hp1 + hp2 - 360) / 2

    T = (
        1
        - 0.17 * numpy.cos((Hpbar - 30) / 180 * numpy.pi)
        + 0.24 * numpy.cos(2 * Hpbar / 180 * numpy.pi)
        + 0.32 * numpy.cos((3 * Hpbar + 6) / 180 * numpy.pi)
        + 0.20 * numpy.cos((4 * Hpbar - 63) / 180 * numpy.pi)
    )
    S_L = 1 + 0.015 * (Lpbar - 50) ** 2 / numpy.sqrt(20 + (Lpbar - 50) ** 2)
    S_C = 1 + 0.045 * Cpbar
    S_H = 1 + 0.015 * Cpbar * T
    R_T = (
        -2
        * numpy.sqrt(Cpbar ** 7 / (Cpbar ** 7 + 25 ** 7))
        * numpy.sin(60 * numpy.exp(-(((Hpbar - 275) / 25) ** 2)) / 180 * numpy.pi)
    )

    dE00 = numpy.sqrt(
        (dLp / k_L / S_L) ** 2
        + (dCp / k_C / S_C) ** 2
        + (dHp / k_H / S_H) ** 2
        + R_T * dCp / k_C / S_C * dHp / k_H / S_H
    )

    return dE00
