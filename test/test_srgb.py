# -*- coding: utf-8 -*-
#
import numpy
import pytest

import colorio


@pytest.mark.parametrize(
    "vals", [numpy.random.rand(3), numpy.random.rand(3, 7), numpy.random.rand(3, 4, 5)]
)
def test_conversion(vals):
    srgb_linear = colorio.SrgbLinear()

    out = srgb_linear.to_xyz100(srgb_linear.from_xyz100(vals))
    assert numpy.all(abs(vals - out) < 1.0e-14)

    out = srgb_linear.to_srgb1(srgb_linear.from_srgb1(vals))
    assert numpy.all(abs(vals - out) < 1.0e-14)
    return


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
    srgb_linear = colorio.SrgbLinear()
    assert numpy.all(
        abs(srgb_linear.from_srgb1(vals) - ref) < 1.0e-14 * numpy.array(ref)
    )
    return


@pytest.mark.parametrize(
    "vals,ref",
    [
        ([0.00650735931, 0.00789021442, 0.114259116060], [2.61219, 1.52732, 10.96471]),
        ([0.03836396959, 0.01531740787, 0.014587362033], [2.39318, 2.01643, 1.64315]),
    ],
)
def test_reference_xyz(vals, ref):
    srgb_linear = colorio.SrgbLinear()
    assert numpy.all(abs(srgb_linear.to_xyz100(vals) - ref) < 1.0e-3 * numpy.array(ref))
    return


def test_whitepoint():
    srgb_linear = colorio.SrgbLinear()
    val = srgb_linear.to_xyz100([1.0, 1.0, 1.0])
    d65_whitepoint = colorio.illuminants.whitepoints_cie1931["D65"]
    assert numpy.all(numpy.abs(val - d65_whitepoint) < 1.0e-12)
    return
