# -*- coding: utf-8 -*-
#
import colorio

import matplotlib.pyplot as plt
import numpy
import pytest


@pytest.mark.parametrize('illuminant,decimals,values', [
    (colorio.illuminants.a(5), 5, [0.93048, 1.12821, 1.35769]),
    (colorio.illuminants.d50(), 3, [0.019, 2.051, 7.778]),
    (colorio.illuminants.d55(), 3, [0.024, 2.072, 11.224]),
    (colorio.illuminants.d65(), 4, [0.03410, 3.2945, 20.2360]),
    # 5.132 is different from the standard; 5.133 is listed there. This is a
    # false rounding.
    (colorio.illuminants.d75(), 3, [0.043, 5.132, 29.808]),
    ])
def test_values(illuminant, decimals, values):
    lmbda, data = illuminant
    rdata = numpy.around(data, decimals=decimals)
    assert rdata[0] == values[0]
    assert rdata[1] == values[1]
    assert rdata[2] == values[2]
    return


def test_white_point():
    values = colorio.illuminants.white_point(colorio.illuminants.d65())
    assert numpy.all(abs(values - [0.9504716, 1.0, 1.08883546]) < 1.0e-8)
    return


def test_show():
    lmbda, data = colorio.illuminants.d65()
    plt.plot(lmbda, data)
    plt.ylim(ymin=0)
    plt.show()
    return


if __name__ == '__main__':
    test_white_point()
    # test_show()
