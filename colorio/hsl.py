import numpy


class Hsl(object):
    def from_srgb1(self, srgb1):
        srgb = numpy.asarray(srgb1, dtype=float)
        assert numpy.all(srgb >= 0)
        assert numpy.all(srgb <= 1)

        argmax = numpy.argmax(srgb)
        max_val = srgb[argmax]
        argmin = numpy.argmin(srgb)
        min_val = srgb[argmin]

        diff = max_val - min_val

        if max_val == min_val:
            H = 0.0
        else:
            if argmax == 0:
                H = 60 * (srgb[1] - srgb[2]) / diff
            elif argmax == 1:
                H = 60 * (2 + ((srgb[2] - srgb[0]) / diff))
            else:
                assert argmax == 2
                H = 60 * (4 + (srgb[0] - srgb[1]) / diff)

            if H < 0:
                H += 360

        if max_val == 0 or min_val == 1:
            S = 0
        else:
            S = diff / (1 - numpy.abs(max_val + min_val - 1))

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
