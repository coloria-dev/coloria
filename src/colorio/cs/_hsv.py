"""
https://en.wikipedia.org/wiki/HSL_and_HSV
"""
import numpy as np


class HSV:
    def from_srgb1(self, srgb1):
        srgb = np.asarray(srgb1, dtype=float)
        orig_shape = srgb.shape
        srgb = srgb.reshape(3, -1)
        assert np.all(srgb >= 0)
        assert np.all(srgb <= 1)

        argmax = np.argmax(srgb, axis=0)
        max_val = np.max(srgb, axis=0)
        min_val = np.min(srgb, axis=0)

        diff = max_val - min_val

        H = np.empty(srgb.shape[1:], dtype=float)
        H[max_val == min_val] = 0
        i = argmax == 0
        H[i] = 60 * (0 + (srgb[1][i] - srgb[2][i]) / diff[i])
        i = argmax == 1
        H[i] = 60 * (2 + (srgb[2][i] - srgb[0][i]) / diff[i])
        i = argmax == 2
        H[i] = 60 * (4 + (srgb[0][i] - srgb[1][i]) / diff[i])
        H = np.mod(H, 360)

        S = np.empty(srgb.shape[1:], dtype=float)
        S[max_val == 0] = 0
        i = (max_val > 0) & (min_val < 1)
        S[i] = diff[i] / max_val[i]

        H = H.reshape(orig_shape[1:])
        S = S.reshape(orig_shape[1:])
        V = max_val.reshape(orig_shape[1:])
        return np.array([H, S, V])

    def to_srgb1(self, hsl):
        hsl = np.asarray(hsl)
        H, S, V = hsl
        if not np.all((0 <= H) & (H <= 360)):
            raise ValueError("Illegal values in H.")
        if not np.all((0 <= S) & (S <= 1)):
            raise ValueError("Illegal values in S.")
        if not np.all((0 <= V) & (V <= 1)):
            raise ValueError("Illegal values in V.")

        C = V * S
        H_dash = H / 60
        X = C * (1 - np.abs(np.mod(H_dash, 2) - 1))
        Z = np.zeros(C.shape)

        R1 = np.empty(C.shape)
        G1 = np.empty(C.shape)
        B1 = np.empty(C.shape)
        i = (0 <= H_dash) & (H_dash <= 1)
        R1[i], G1[i], B1[i] = C[i], X[i], Z[i]
        i = (1 < H_dash) & (H_dash <= 2)
        R1[i], G1[i], B1[i] = X[i], C[i], Z[i]
        i = (2 < H_dash) & (H_dash <= 3)
        R1[i], G1[i], B1[i] = Z[i], C[i], X[i]
        i = (3 < H_dash) & (H_dash <= 4)
        R1[i], G1[i], B1[i] = Z[i], X[i], C[i]
        i = (4 < H_dash) & (H_dash <= 5)
        R1[i], G1[i], B1[i] = X[i], Z[i], C[i]
        i = (5 < H_dash) & (H_dash <= 6)
        R1[i], G1[i], B1[i] = C[i], Z[i], X[i]

        m = V - C
        return np.array([R1 + m, G1 + m, B1 + m])

    def from_srgb256(self, srgb256):
        return self.from_srgb1(np.asarray(srgb256) / 255.0)
