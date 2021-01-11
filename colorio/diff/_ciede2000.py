import numpy


def ciede2000(lab1, lab2, k_L=1.0, k_C=1.0, k_H=1.0):
    lab1 = numpy.asarray(lab1)
    lab2 = numpy.asarray(lab2)

    L1, a1, b1 = lab1
    L2, a2, b2 = lab2

    C1 = numpy.sqrt(a1 ** 2 + b1 ** 2)
    C2 = numpy.sqrt(a2 ** 2 + b2 ** 2)
    C_mean = (C1 + C2) / 2

    G = 0.5 * (1 - numpy.sqrt(C_mean ** 7 / (C_mean ** 7 + 25 ** 7)))
    a1p = (1 + G) * a1
    a2p = (1 + G) * a2

    C1p = numpy.sqrt(a1p ** 2 + b1 ** 2)
    C2p = numpy.sqrt(a2p ** 2 + b2 ** 2)

    h1p = numpy.degrees(numpy.arctan2(b1, a1p)) % 360
    h2p = numpy.degrees(numpy.arctan2(b2, a2p)) % 360
    # make sure dhp is in the (-180, 180) range
    hp_diff = h2p - h1p
    dhp = ((hp_diff + 180) % 360) - 180

    dLp = L2 - L1
    dCp = C2p - C1p
    dHp = 2 * numpy.sqrt(C1p * C2p) * numpy.sin(numpy.radians(dhp / 2))

    Lp_mean = (L1 + L2) / 2
    Cp_mean = (C1p + C2p) / 2

    # TODO clean this up
    hp_mean = numpy.empty_like(h1p)
    hp_avg = (h1p + h2p) / 2
    idx = numpy.abs(hp_diff) <= 180
    hp_mean[idx] = hp_avg[idx]
    idx = ~(numpy.abs(hp_diff) <= 180) & (hp_avg < 180)
    hp_mean[idx] = hp_avg[idx] + 180
    idx = ~(numpy.abs(hp_diff) <= 180) & (hp_avg >= 180)
    hp_mean[idx] = hp_avg[idx] - 180
    idx = (C1p == 0.0) | (C2p == 0.0)
    hp_mean[idx] = 2 * hp_avg[idx]

    T = (
        1.0
        - 0.17 * numpy.cos(numpy.radians(hp_mean - 30))
        + 0.24 * numpy.cos(numpy.radians(2 * hp_mean))
        + 0.32 * numpy.cos(numpy.radians(3 * hp_mean + 6))
        - 0.20 * numpy.cos(numpy.radians(4 * hp_mean - 63))
    )
    dtheta = 30 * numpy.exp(-(((hp_mean - 275) / 25) ** 2))
    R_C = 2 * numpy.sqrt(Cp_mean ** 7 / (Cp_mean ** 7 + 25 ** 7))
    S_L = 1 + 0.015 * (Lp_mean - 50) ** 2 / numpy.sqrt(20 + (Lp_mean - 50) ** 2)
    S_C = 1 + 0.045 * Cp_mean
    S_H = 1 + 0.015 * Cp_mean * T
    R_T = -numpy.sin(numpy.radians(2 * dtheta)) * R_C

    dE00 = numpy.sqrt(
        (dLp / k_L / S_L) ** 2
        + (dCp / k_C / S_C) ** 2
        + (dHp / k_H / S_H) ** 2
        + R_T * (dCp / k_C / S_C) * (dHp / k_H / S_H)
    )

    return dE00
