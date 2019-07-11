import numpy


class Hsl(object):
    def from_srgb1(self, srgb1):
        srgb = numpy.asarray(srgb1, dtype=float)
        assert numpy.all(srgb >= 0)
        assert numpy.all(srgb <= 1)

        argmax = numpy.argmax(srgb, axis=0)
        max_val = numpy.max(srgb, axis=0)
        min_val = numpy.min(srgb, axis=0)

        diff = max_val - min_val

        H = numpy.empty(srgb1.shape[1:], dtype=float)
        H[max_val == min_val] = 0
        i = argmax == 0
        H[i] = 60 * (0 + (srgb[1][i] - srgb[2][i]) / diff[i])
        i = argmax == 1
        H[i] = 60 * (2 + (srgb[2][i] - srgb[0][i]) / diff[i])
        i = argmax == 2
        H[i] = 60 * (4 + (srgb[0][i] - srgb[1][i]) / diff[i])
        H = numpy.mod(H, 360)

        S = numpy.empty(srgb1.shape[1:], dtype=float)
        S[max_val == 0] = 0
        S[min_val == 1] = 0
        S[(max_val > 0) & (min_val < 1)] = diff / (1 - numpy.abs(max_val + min_val - 1))

        L = (max_val + min_val) / 2

        return numpy.array([H, S, L])

    def to_srgb1(self, hsl):
        H, S, L = hsl
        assert numpy.all(H >= 0)
        assert numpy.all(H <= 360)
        assert numpy.all(S >= 0)
        assert numpy.all(S <= 1)
        assert numpy.all(L >= 0)
        assert numpy.all(L <= 1)

        C = (1 - numpy.abs(2 * L - 1)) * S
        H_dash = H / 60
        X = C * (1 - numpy.abs(numpy.mod(H_dash, 2) - 1))
        if H_dash <= 1:
            R1, G1, B1 = C, X, 0
        elif H_dash <= 2:
            R1, G1, B1 = X, C, 0
        elif H_dash <= 3:
            R1, G1, B1 = 0, C, X
        elif H_dash <= 4:
            R1, G1, B1 = 0, X, C
        elif H_dash <= 5:
            R1, G1, B1 = X, 0, C
        else:
            assert H_dash <= 6
            R1, G1, B1 = C, 0, X

        m = L - C / 2
        return numpy.array([R1 + m, G1 + m, B1 + m])

    def from_srgb256(self, srgb256):
        return self.from_srgb1(numpy.asarray(srgb256) / 255.0)
