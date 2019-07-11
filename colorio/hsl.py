import numpy


class Hsl(object):
    def from_srgb1(self, srgb1):
        srgb = numpy.asarray(srgb1, dtype=float)

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
            S = diff / (1 - numpy.abs(max_val + min_val - 1)) * 100

        L = (max_val + min_val) / 2 * 100

        return numpy.array([H, S, L])

    def to_srgb1(self, srgb_linear):
        a = 0.055
        is_smaller = srgb_linear <= 0.0031308
        srgb = numpy.array(srgb_linear, dtype=float)
        srgb[is_smaller] *= 12.92
        srgb[~is_smaller] = (1 + a) * srgb[~is_smaller] ** (1 / 2.4) - a
        return srgb

    def from_srgb256(self, srgb256):
        return self.from_srgb1(numpy.asarray(srgb256) / 255.0)
