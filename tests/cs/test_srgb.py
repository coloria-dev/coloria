import numpy as np
import pytest

import colorio
from colorio.cs import ColorCoordinates, Srgb1, Srgb255, SrgbLinear, convert

rng = np.random.default_rng(0)


@pytest.mark.parametrize("cs", [SrgbLinear(), Srgb1(), Srgb255()])
@pytest.mark.parametrize(
    "vals", [100 * rng.random(3), 100 * rng.random((3, 7)), 100 * rng.random((3, 4, 5))]
)
def test_conversion(cs, vals):
    out = cs.to_xyz100(cs.from_xyz100(vals, mode="ignore"))
    assert np.all(abs(vals - out) < 1.0e-13)


@pytest.mark.parametrize(
    "vals,ref",
    [
        ([0.1, 0.2, 0.3], [0.010022825574869, 0.0331047665708851, 0.0732389558784054]),
        ([0.9, 0.8, 0.7], [0.787412289395617, 0.603827338855338, 0.447988412441883]),
        (
            [0.04, 0.02, 0.01],
            [0.00309597523219814, 0.00154798761609907, 0.000773993808049536],
        ),
    ],
)
def test_reference_srgb(vals, ref):
    vals = ColorCoordinates(vals, Srgb1())
    assert np.all(abs(convert(vals, SrgbLinear()).data - ref) < 1.0e-14 * np.array(ref))


@pytest.mark.parametrize(
    "vals,ref",
    [
        ([0.00650735931, 0.00789021442, 0.114259116060], [2.61219, 1.52732, 10.96471]),
        ([0.03836396959, 0.01531740787, 0.014587362033], [2.39318, 2.01643, 1.64315]),
    ],
)
def test_reference_xyz(vals, ref):
    srgb_linear = SrgbLinear()
    assert np.all(abs(srgb_linear.to_xyz100(vals) - ref) < 1.0e-3 * np.array(ref))


def test_whitepoint():
    srgb_linear = SrgbLinear()
    val = srgb_linear.to_xyz100([1.0, 1.0, 1.0])
    d65_whitepoint = colorio.illuminants.whitepoints_cie1931["D65"]
    assert np.all(np.abs(val - d65_whitepoint) < 1.0e-12)


def test_modes():
    xyz = [83.0, 53.0, 67.0]
    srgb_linear = colorio.cs.SrgbLinear()

    with pytest.raises(Exception):
        srgb_linear.from_xyz100(xyz, mode="error")

    rgb = srgb_linear.from_xyz100(xyz, mode="nan")
    assert np.isnan(rgb[0])

    rgb = srgb_linear.from_xyz100(xyz, mode="ignore")
    ref = [1.5411487959020491, 0.21767779000754928, 0.6463796478923135]
    assert np.all(np.abs(rgb - ref) < 1.0e-13 * np.abs(ref))

    rgb = srgb_linear.from_xyz100(xyz, mode="clip")
    ref = [1.0, 0.21767779000754928, 0.6463796478923135]
    assert np.all(np.abs(rgb - ref) < 1.0e-13 * np.abs(ref))
