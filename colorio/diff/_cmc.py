import numpy as np

from ..cs import CIELAB, CIELCH


# Acceptability: l=2, c=1 Perceptability: l=1, c=1
def cmc(lab1, lab2, l=2.0, c=1.0):
    cielab = CIELAB()
    cielch = CIELCH()
    lch1 = cielch.from_xyz100(cielab.to_xyz100(lab1))
    lch2 = cielch.from_xyz100(cielab.to_xyz100(lab2))

    L1, C1, h1 = lch1
    L2, C2, _ = lch2

    F = np.sqrt(C1 ** 4 / (C1 ** 4 + 1900))
    if 164 <= h1 and h1 <= 345:
        T = 0.56 + np.abs(0.2 * np.cos(np.radians(h1 + 168)))
    else:
        T = 0.36 + np.abs(0.4 * np.cos(np.radians(h1 + 35)))

    if L1 < 16:
        S_L = 0.511
    else:
        S_L = (0.040975 * L1) / (1 + 0.01765 * L1)

    S_C = 0.0638 * C1 / (1 + 0.0131 * C1) + 0.638
    S_H = S_C * (F * T + 1 - F)

    L1, a1, b1 = lab1
    L2, a2, b2 = lab2
    dC = C1 - C2
    da = a1 - a2
    db = b1 - b2
    dHab2 = da ** 2 + db ** 2 - dC ** 2

    dE = np.sqrt(
        ((L2 - L1) / l / S_L) ** 2 + (dC / c / S_C) ** 2 + dHab2 / S_H ** 2
    )
    return dE
