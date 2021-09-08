import numpy as np
import pytest

import colorio
from colorio.illuminants import whitepoints_cie1931

ob2 = colorio.observers.cie_1931_2()


@pytest.mark.parametrize(
    "illuminant,decimals,values",
    [
        (colorio.illuminants.a(5), 5, [0.93048, 1.35769, 1.92508]),
        (colorio.illuminants.d50(), 3, [0.019, 2.051, 7.778]),
        (colorio.illuminants.d55(), 3, [0.024, 2.072, 11.224]),
        (colorio.illuminants.d65(), 4, [0.03410, 3.2945, 20.2360]),
        # 5.132 is different from the standard; 5.133 is listed there.
        # This is a false rounding.
        (colorio.illuminants.d75(), 3, [0.043, 5.132, 29.808]),
    ],
)
def test_values(illuminant, decimals, values):
    rdata = np.around(illuminant.data, decimals=decimals)
    print(illuminant)
    print(list(rdata[:6:2]))
    assert np.all(rdata[:6:2] == values)


# Don't try to match compute_whitepoint() with the actual whitepoints given in the
# standard. Too many things are not clear:
#
#   - What's the underlying data?
#   - How are the values interpolated?
#   - How are the products integrated?
#
# The methods for interpolation/integration given in the standard are quite dubious,
# too. Don't get into that.
#
@pytest.mark.parametrize(
    "illuminant,observer,ref",
    [
        (colorio.illuminants.a(), ob2, whitepoints_cie1931["A"]),
        (colorio.illuminants.d50(), ob2, [96.424, 100.0, 82.513]),
        (colorio.illuminants.d65(), ob2, [95.047, 100.0, 108.883]),
        (colorio.illuminants.d75(), ob2, [94.972, 100.0, 122.619]),
        # Should all sum up to 100, but the standard doesn't actually provide values
        # that do so. There's a newer not-yet-standard, CIE 2012, where this condition
        # is fulfilled up to approx 1.0e-6.
        (colorio.illuminants.e(), ob2, [100.008, 100.0, 100.033]),
        (colorio.illuminants.f2(), ob2, [103.755, 100.0, 49.864]),
        (colorio.illuminants.f7(), ob2, [95.019, 100.0, 108.639]),
        (colorio.illuminants.f11(), ob2, [100.904, 100.0, 64.284]),
    ],
)
def test_white_point(illuminant, observer, ref):
    print(illuminant)
    print(observer)
    values = colorio.illuminants.compute_whitepoint(illuminant, observer)
    print(list(values))
    print(list(np.round(values, 3)))
    assert np.all(np.round(values, 3) == ref)


def test_spectrum_to_xyz100():
    spectrum = colorio.illuminants.d65()
    observer = colorio.observers.cie_1931_2()
    out = colorio.illuminants.spectrum_to_xyz100(spectrum, observer)
    out = out / out[1] * 100
    print(list(out))
    ref = [95.04698658153113, 100.0, 108.88271743218647]
    assert np.all(np.abs(out - ref) < 1.0e-13 * np.abs(ref))


if __name__ == "__main__":
    test_white_point(
        colorio.illuminants.d50(), ob2, colorio.illuminants.whitepoints_cie1931["D50"]
    )
