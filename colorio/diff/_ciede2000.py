import numpy


def ciede2000(lab1, lab2, k_L=1.0, k_C=1.0, k_H=1.0):
    L1, a1, b1 = lab1
    L2, a2, b2 = lab2

    C1 = numpy.sqrt(a1 ** 2 + b1 ** 2)
    C2 = numpy.sqrt(a2 ** 2 + b2 ** 2)
    Cbar = (C1 + C2) / 2

    G = 0.5 * (1 - numpy.sqrt(Cbar ** 7 / (Cbar ** 7 + 25 ** 7)))
    a1p = (1 + G) * a1
    a2p = (1 + G) * a2

    C1p = numpy.sqrt(a1p ** 2 + b1 ** 2)
    C2p = numpy.sqrt(a2p ** 2 + b2 ** 2)

    h1p = numpy.mod(numpy.degrees(numpy.arctan2(b1, a1p)), 360)
    h2p = numpy.mod(numpy.degrees(numpy.arctan2(b2, a2p)), 360)

    dLp = L2 - L1
    dCp = C2p - C1p

    if C1p == 0.0 or C2p == 0.0:
        dhp = 0.0
    else:
        if abs(h2p - h1p) <= 180:
            dhp = h2p - h1p
        elif h2p - h1p > 180:
            dhp = h2p - h1p - 360
        else:
            assert h2p - h1p < -180
            dhp = h2p - h1p + 360

    dHp = 2 * numpy.sqrt(C1p * C2p) * numpy.sin(numpy.radians(dhp / 2))

    Lpbar = (L1 + L2) / 2
    Cpbar = (C1p + C2p) / 2

    if C1p == 0.0 or C2p == 0.0:
        Hpbar = h1p + h2p
    elif numpy.abs(h1p - h2p) <= 180:
        Hpbar = (h1p + h2p) / 2
    else:
        if h1p + h2p < 360:
            Hpbar = (h1p + h2p + 360) / 2
        else:
            Hpbar = (h1p + h2p - 360) / 2

    T = (
        1.0
        - 0.17 * numpy.cos(numpy.radians(Hpbar - 30))
        + 0.24 * numpy.cos(numpy.radians(2 * Hpbar))
        + 0.32 * numpy.cos(numpy.radians(3 * Hpbar + 6))
        - 0.20 * numpy.cos(numpy.radians(4 * Hpbar - 63))
    )
    dtheta = 30 * numpy.exp(-((Hpbar - 275) / 25) ** 2)
    R_C = 2 * numpy.sqrt(Cpbar ** 7 / (Cpbar ** 7 + 25 ** 7))
    S_L = 1 + 0.015 * (Lpbar - 50) ** 2 / numpy.sqrt(20 + (Lpbar - 50) ** 2)
    S_C = 1 + 0.045 * Cpbar
    S_H = 1 + 0.015 * Cpbar * T
    R_T = -numpy.sin(numpy.radians(2 * dtheta)) * R_C

    dE00 = numpy.sqrt(
        (dLp / k_L / S_L) ** 2
        + (dCp / k_C / S_C) ** 2
        + (dHp / k_H / S_H) ** 2
        + R_T * (dCp / k_C / S_C) * (dHp / k_H / S_H)
    )

    return dE00
