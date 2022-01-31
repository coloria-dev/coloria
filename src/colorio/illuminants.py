import json
import pathlib

import numpy as np

from ._helpers import SpectralData

# The "standard" 2 degree observer (CIE 1931). From
# <https://github.com/njsmith/colorspacious/blob/master/colorspacious/illuminants.py>
whitepoints_cie1931 = {
    "A": np.array([109.850, 100, 35.585]),
    "C": np.array([98.074, 100, 118.232]),
    "D50": np.array([96.422, 100, 82.521]),
    "D55": np.array([95.682, 100, 92.149]),
    "D65": np.array([95.047, 100, 108.883]),
    "D75": np.array([94.972, 100, 122.638]),
    "F2": np.array([99.186, 100, 67.393]),
    "F7": np.array([95.041, 100, 108.747]),
    "F11": np.array([100.962, 100, 64.350]),
}

# The "supplementary" 10 degree observer (CIE 1964). From
# <https://github.com/njsmith/colorspacious/blob/master/colorspacious/illuminants.py>
whitepoints_cie1964 = {
    "A": np.array([111.144, 100, 35.200]),
    "C": np.array([97.285, 100, 116.145]),
    "D50": np.array([96.720, 100, 81.427]),
    "D55": np.array([95.799, 100, 90.926]),
    "D65": np.array([94.811, 100, 107.304]),
    "D75": np.array([94.416, 100, 120.641]),
    "F2": np.array([103.279, 100, 69.027]),
    "F7": np.array([95.792, 100, 107.686]),
    "F11": np.array([103.863, 100, 65.607]),
}


def spectrum_to_xyz100(
    spectrum: SpectralData,
    observer: SpectralData,
    interpolation_type: str = "linear",
) -> np.ndarray:
    """Computes the tristimulus values XYZ from a given spectrum for a given observer
    via

      X_i = sum_lambda spectrum_i(lambda) * observer_i(lambda) delta_lambda.

    The technical report CIE Standard Illuminants for Colorimetry, 1999, section 7
    ("Recommendations concerning the calculation of tristimulus values and chromaticity
    coordinates"), gives a recommendation on how to perform the computation:

    > The CIE Standard (CIE, 1986a) on standard colorimetric observers recommends that
    > the CIE tristimulus values of a colour stimulus be obtained by multiplying at each
    > wavelength the value of the colour stimulus function phi_lambda(lambda) by that of
    > each of the CIE colour-matching functions and integrating each set of products
    > over the wavelength range corresponding to the entire visible spectrum, 360 nm to
    > 830 nm. The integration can be carried out by numerical summation at wavelength
    > intervals, delta lmbda, equal to 1 nm.

    Note that the above sum is supposed to approximate the integral

      int_lambda spectrum_i(lambda) * observer_i(lambda) dlambda.

    It gets a little tricky when asking what the correct value is for monochromatic
    light. Mathematically, one would insert a scaled delta distribution for spectrum_i
    to get a * observer_i(lambda0). It would be easy to assume that you do the same in
    the corresponding sum above, but there you get

      a * observer_i(lambda0) * delta_lambda.

    So just using a unity vector for spectrum_i in the sum does _not_ correspond to
    monochromatic light, but rather light of a very small bandwidth (between
    lambda0 - delta_lambda and lambda0 + delta_lambda).

    To get a more consistent view on things, it is useful to look at observer_i as a
    piecewise linear function, and the spectrum vector as a piecewise constant function.
    The above sum accurately represents the integral of the product of the two.

    Note that any constant factor (like delta_lambda) gets canceled out in x and y (of
    xyY), so being careless might not be punished in all applications.
    """
    # The technical document prescribes that the integration be performed over
    # the wavelength range corresponding to the entire visible spectrum, 360 nm
    # to 830 nm. Make sure the observer has the appropriate data.
    delta = 1
    lmbda = np.arange(360, 831, delta)
    assert np.all(observer.lmbda_nm == lmbda)

    # Adapt the spectrum
    mask = (360 <= spectrum.lmbda_nm) & (spectrum.lmbda_nm <= 830)
    lambda_s = spectrum.lmbda_nm[mask]
    data_s = spectrum.data[mask]

    if not np.array_equal(lambda_s, lmbda):
        # The technical report specifies the interpolation techniques, too:
        # ```
        # Use one of the four following methods to calculate needed but unmeasured
        # values of phi(l), R(l) or tau(l) within the range of measurements:
        #   1) the third-order polynomial interpolation (Lagrange) from the four
        #      neighbouring data points around the point to be interpolated, or
        #   2) cubic spline interpolation formula, or
        #   3) a fifth order polynomial interpolation formula from the six
        #      neighboring data points around the point to be interpolated, or
        #   4) a Sprague interpolation (see Seve, 2003).
        # ```
        if interpolation_type == "linear":
            # Linear interpolation isn't actually part of the standard. All types
            # except this one are bad. One reason: Results can be negative even if the
            # function isn't.
            data_s = np.interp(observer.lmbda_nm, lambda_s, data_s)
        if interpolation_type == "lagrange-3":
            import scipyx

            poly = scipyx.interp_rolling_lagrange(lambda_s, data_s, order=3)
            data_s = poly(lmbda)
        elif interpolation_type == "cubic spline":
            # The standard doesn't give the boundary conditions
            from scipy.interpolate import CubicSpline

            cs = CubicSpline(lambda_s, data_s, bc_type="not-a-knot")
            data_s = cs(lmbda)
        elif interpolation_type == "lagrange-5":
            assert interpolation_type == "lagrange-5"
            import scipyx

            poly = scipyx.interp_rolling_lagrange(lambda_s, data_s, order=5)
            data_s = poly(lmbda)

    xyz100 = np.sum(data_s * observer.data * delta, axis=1)

    # For transmittant or reflectant objects, data_s is typically the illuminant
    # spectrum S times a reflectance/transmittance factor R, 0<=R<=1. In this case, one
    # should multiply by
    #
    # k = 100 / np.sum(S * observer.data[1] * delta)
    #

    return xyz100


def compute_whitepoint(illuminant: SpectralData, observer: SpectralData) -> np.ndarray:
    xyz100 = spectrum_to_xyz100(illuminant, observer)
    # make sure the Y value is 100
    xyz100 *= 100 / xyz100[1]
    return xyz100


def planckian_radiator(temperature):
    lmbda_nm = np.arange(300, 831)
    # light speed
    c = 299792458.0
    # Plank constant
    h = 6.62607015e-34
    # Boltzmann constant
    k = 1.380649e-23
    c1 = 2 * np.pi * h * c**2
    c2 = h * c / k
    lmbda = 1.0e-9 * lmbda_nm
    return SpectralData(
        lmbda_nm,
        c1 / lmbda**5 / (np.exp(c2 / lmbda / temperature) - 1),
        f"Planckian radiator ({temperature} K)",
    )


def a(interval_nm: int = 1):
    """CIE Standard Illuminants for Colorimetry, 1999:
    CIE standard illuminant A is intended to represent typical, domestic,
    tungsten-filament lighting. Its relative spectral power distribution is that of a
    Planckian radiator at a temperature of approximately 2856 K. CIE standard illuminant
    A should be used in all applications of colorimetry involving the use of
    incandescent lighting, unless there are specific reasons for using a different
    illuminant.
    """
    # https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_A
    lmbda_nm = np.arange(300, 831, interval_nm)
    # When Illuminant A was standardized, the natural constants where such that this was
    # the value of c2. The values of the constants have since been revised. In order to
    # avoid further possible changes in the color temperature, the CIE now specifies the
    # SPD directly, based on the original (1931) value of c2.
    c2 = 1.435e-2
    color_temp = 2848
    vals = (
        100
        * (560 / lmbda_nm) ** 5
        * (
            (np.exp(c2 / (color_temp * 560e-9)) - 1)
            / (np.exp(c2 / (color_temp * lmbda_nm * 1.0e-9)) - 1)
        )
    )
    return SpectralData(lmbda_nm, vals, "Illuminant A")


def c():
    this_dir = pathlib.Path(__file__).resolve().parent
    return _from_file(this_dir / "data/illuminants/c.json")


def d(nominal_temperature: float):
    """CIE D-series illuminants.

    The technical report `Colorimetry, 3rd edition, 2004` gives the data for D50, D55,
    and D65 explicitly, but also explains how it's computed for S0, S1, S2. Values are
    given at 5nm resolution in the document, but really every other value is just
    interpolated. Hence, only provide 10 nm data here.
    """
    # From CIE 15:2004. Colorimetry, 3rd edition, 2004 (page 69, note 5):
    #
    # The method required to calculate the values for the relative spectral power
    # distributions of illuminants D50, D55, D65, and D75, in Table T.1 is as follows
    #   1. Multiply the nominal correlated colour temperature (5000 K, 5500 K, 6500 K or
    #   7500 K) by 1,4388/1,4380.
    #   2. Calculate XD and YD using the equations given in the text.
    #   3. Calculate M1 and M2 using the equations given in the text.
    #   4. Round M1 and M2 to three decimal places.
    #   5. Calculate S(lambda) every 10 nm by
    #        S(lambda) = S0(lambda) + M1 S1(lambda) + M2 S2(lambda)
    #      using values of S0(lambda), S1(lambda) and S2(lambda) from Table T.2.
    #   6. Interpolate the 10 nm values of S(lambda) linearly to obtain values at
    #      intermediate wavelengths.
    tcp = 1.4388e-2 / 1.4380e-2 * nominal_temperature

    if 4000 <= tcp <= 7000:
        xd = ((-4.6070e9 / tcp + 2.9678e6) / tcp + 0.09911e3) / tcp + 0.244063
    else:
        assert 7000 < tcp <= 25000
        xd = ((-2.0064e9 / tcp + 1.9018e6) / tcp + 0.24748e3) / tcp + 0.237040

    yd = (-3.000 * xd + 2.870) * xd - 0.275

    m1 = (-1.3515 - 1.7703 * xd + 5.9114 * yd) / (0.0241 + 0.2562 * xd - 0.7341 * yd)
    m2 = (+0.0300 - 31.4424 * xd + 30.0717 * yd) / (0.0241 + 0.2562 * xd - 0.7341 * yd)

    m1 = np.around(m1, decimals=3)
    m2 = np.around(m2, decimals=3)

    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "data/illuminants/d.json") as f:
        data = json.load(f)

    # https://www.wikiwand.com/en/Standard_illuminant:
    # > The tabulated SPDs presented by the CIE today are derived by linear
    # > interpolation of the 10 nm data set down to 5 nm.
    lmbda_start, lmbda_end, lmbda_step = data["lambda_nm"]
    assert lmbda_step == 10
    lmbda10 = np.arange(lmbda_start, lmbda_end + 1, lmbda_step)
    S10 = np.asarray(data["S"])
    # 6. Interpolate the 10 nm values of S(lambda) linearly to obtain values at
    #    intermediate wavelengths.
    # nschloe: I think that's pretty useless, but yeah, that's the standard.
    lmbda5 = np.arange(lmbda_start, lmbda_end + 1, 5)
    S5 = np.array(
        [
            np.interp(lmbda5, lmbda10, S10[0]),
            np.interp(lmbda5, lmbda10, S10[1]),
            np.interp(lmbda5, lmbda10, S10[2]),
        ]
    )

    return SpectralData(
        lmbda5,
        S5[0] + m1 * S5[1] + m2 * S5[2],
        "Illuminant D" + str(nominal_temperature)[:2],
    )


def d50():
    """CIE illuminant D50, mid-morning/mid-afternoon daylight, at 10nm resolution."""
    return d(5000)


def d55():
    """CIE illuminant D55, mid-morning/mid-afternoon daylight, at 10nm resolution."""
    return d(5500)


def d65():
    """CIE standard illuminant D65, sampled at 10nm intervals."""
    return d(6500)


def d75():
    """CIE illuminant D75"""
    return d(7500)


def e():
    """This is a hypothetical reference radiator. All wavelengths in CIE illuminant E
    are weighted equally with a relative spectral power of 100.0.
    """
    return SpectralData(np.arange(300, 831), np.full(531, 100.0), "Illuminant E")


def _from_file(filename):
    with open(filename) as f:
        data = json.load(f)

    start, stop, step = data["lambda_nm"]
    return SpectralData(
        np.arange(start, stop + 1, step),
        data["values"],
        data["description"],
    )


def f2():
    this_dir = pathlib.Path(__file__).resolve().parent
    return _from_file(this_dir / "data/illuminants/f2.json")


def f7():
    this_dir = pathlib.Path(__file__).resolve().parent
    return _from_file(this_dir / "data/illuminants/f7.json")


def f11():
    this_dir = pathlib.Path(__file__).resolve().parent
    return _from_file(this_dir / "data/illuminants/f11.json")
