"""
https://en.wikipedia.org/wiki/HSL_and_HSV
"""
import numpy as np


class HSL:
    def from_rgb1(self, srgb1):
        srgb = np.asarray(srgb1, dtype=float)
        orig_shape = srgb.shape
        srgb = srgb.reshape(3, -1)
        assert not np.any(srgb < 0.0)
        assert not np.any(srgb > 1.0)

        argmax = np.argmax(srgb, axis=0)
        max_val = np.max(srgb, axis=0)
        min_val = np.min(srgb, axis=0)

        diff = max_val - min_val

        is_diff_0 = diff == 0.0

        H = np.empty(srgb.shape[1:], dtype=float)
        # Set hue to 0.0 for grey values. Could do anything here, nan would be
        # reasonable too.
        H[is_diff_0] = 0.0
        i = (argmax == 0) & ~is_diff_0
        H[i] = 60 * (0 + (srgb[1][i] - srgb[2][i]) / diff[i])
        i = (argmax == 1) & ~is_diff_0
        H[i] = 60 * (2 + (srgb[2][i] - srgb[0][i]) / diff[i])
        i = (argmax == 2) & ~is_diff_0
        H[i] = 60 * (4 + (srgb[0][i] - srgb[1][i]) / diff[i])
        H = np.mod(H, 360)

        S = np.empty(srgb.shape[1:], dtype=float)
        idx = (max_val > 0) & (min_val < 1)
        S[idx] = diff[idx] / (1 - np.abs(max_val[idx] + min_val[idx] - 1))
        S[~idx] = 0.0 * diff[~idx]

        L = (max_val + min_val) / 2

        H = H.reshape(orig_shape[1:])
        S = S.reshape(orig_shape[1:])
        L = L.reshape(orig_shape[1:])
        return np.array([H, S, L])

    def to_rgb1(self, hsl):
        H, S, L = hsl
        assert not np.any(H < 0)
        assert not np.any(H > 360)
        assert not np.any(S < 0)
        assert not np.any(S > 1)
        assert not np.any(L < 0)
        assert not np.any(L > 1)

        C = (1 - np.abs(2 * L - 1)) * S
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

        m = L - C / 2
        return np.array([R1 + m, G1 + m, B1 + m])

    def from_rgb256(self, srgb256):
        return self.from_rgb1(np.asarray(srgb256) / 255.0)
