# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy
import pytest

import colorio


@pytest.mark.parametrize('illuminant,decimals,values', [
    (colorio.illuminants.a(5e-9), 5, [0.93048, 1.12821, 1.35769]),
    (colorio.illuminants.d50(), 3, [0.019, 2.051, 7.778]),
    (colorio.illuminants.d55(), 3, [0.024, 2.072, 11.224]),
    (colorio.illuminants.d65(), 4, [0.03410, 3.2945, 20.2360]),
    # 5.132 is different from the standard; 5.133 is listed there. This is a
    # false rounding.
    (colorio.illuminants.d75(), 3, [0.043, 5.132, 29.808]),
    ])
def test_values(illuminant, decimals, values):
    _, data = illuminant
    rdata = numpy.around(data, decimals=decimals)
    assert rdata[0] == values[0]
    assert rdata[1] == values[1]
    assert rdata[2] == values[2]
    return


@pytest.mark.parametrize('illuminant,ref,tol', [
    (colorio.illuminants.d65(), [0.95048974, 1.0, 1.08892197], 1.0e-8),
    (colorio.illuminants.e(), [1.00015018, 1.0, 1.00066598], 1.0e-8),
    (colorio.illuminants.f2(), [0.99146841, 1.0, 0.67318498], 1.0e-8),
    ])
def test_white_point(illuminant, ref, tol):
    values = colorio.illuminants.white_point(illuminant)
    print(values)
    print(ref)
    assert numpy.all(abs(values - ref) < tol)
    return


def test_show():
    # lmbda, data = colorio.illuminants.d65()
    for T in [1000, 2000, 3000, 4000, 5000, 1000]:
        lmbda, data = colorio.illuminants.planckian_radiator(T)
        plt.plot(lmbda, data)
    plt.ylim(ymin=0)
    plt.show()
    return


if __name__ == '__main__':
    # test_white_point()
    test_show()
