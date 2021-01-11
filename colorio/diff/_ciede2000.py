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

    h1p = numpy.mod(numpy.degrees(numpy.arctan2(b1, a1p)), 360)
    h2p = numpy.mod(numpy.degrees(numpy.arctan2(b2, a2p)), 360)
    hp_diff = h2p - h1p

    dLp = L2 - L1
    dCp = C2p - C1p

    dhp = numpy.empty_like(h1p)
    idx = numpy.abs(hp_diff) <= 180
    dhp[idx] = hp_diff[idx]
    idx = hp_diff > 180
    dhp[idx] = hp_diff[idx] - 360
    idx = hp_diff < -180
    dhp[idx] = hp_diff[idx] + 360
    #
    is_c_zero = (C1p == 0.0) | (C2p == 0.0)
    dhp[is_c_zero] = 0.0

    dHp = 2 * numpy.sqrt(C1p * C2p) * numpy.sin(numpy.radians(dhp / 2))

    Lp_mean = (L1 + L2) / 2
    Cp_mean = (C1p + C2p) / 2

    hp_mean = numpy.empty_like(h1p)
    hp_sum = h1p + h2p
    idx = is_c_zero
    hp_mean[idx] = hp_sum[idx]
    idx = ~is_c_zero & (numpy.abs(h1p - h2p) <= 180)
    hp_mean[idx] = hp_sum[idx] / 2
    idx = ~is_c_zero & ~(numpy.abs(h1p - h2p) <= 180) & (h1p + h2p < 360)
    hp_mean[idx] = (hp_sum[idx] + 360) / 2
    idx = ~is_c_zero & ~(numpy.abs(h1p - h2p) <= 180) & (h1p + h2p >= 360)
    hp_mean[idx] = (hp_sum[idx] - 360) / 2

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
