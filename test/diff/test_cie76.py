import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "lab1, lab2, ref",
    [
        ((50.0000, 2.6772, -79.7751), (50.0000, 0.0000, -82.7485), 4.001063283678486),
        ((50.0000, -1.1848, -84.8006), (50.0000, 0.0000, -82.7485), 2.3695707311662977),
        ((50.0000, -1.0000, 2.0000), (50.0000, 0.0000, 0.0000), 2.23606797749979),
        ((50.0000, 2.4900, -0.0010), (50.0000, -2.4900, 0.0011), 4.980000442771065),
        ((50.0000, 2.5000, 0.0000), (56.0000, -27.0000, -3.0000), 30.25309901481169),
        ((50.0000, 2.5000, 0.0000), (58.0000, 24.0000, 15.0000), 27.408940147331492),
        ((50.0000, 2.5000, 0.0000), (50.0000, 3.2592, 0.3350), 0.8298250659024465),
        ((61.2901, 3.7196, -5.3901), (61.4292, 2.2480, -4.9620), 1.5389038241553625),
        ((90.9257, -0.5406, -0.9208), (88.6381, -0.8985, -0.7239), 2.3237847964043605),
        ((6.7747, -0.2908, -2.4247), (5.8714, -0.0985, -2.2286), 0.9441320829206047),
        ((2.0776, 0.0795, -1.1350), (0.9033, -0.0636, -0.5514), 1.319108433753647),
    ],
)
def test(lab1, lab2, ref):
    print(lab1, lab2)
    print(ref)
    val = colorio.diff.cie76(lab1, lab2)
    print(val)
    tol = 1.0e-13
    assert np.abs(val - ref) < tol * np.abs(ref)

    # from colormath.color_objects import LabColor
    # from colormath.color_diff import delta_e_cie1976
    # color1 = LabColor(lab_l=lab1[0], lab_a=lab1[1], lab_b=lab1[2])
    # color2 = LabColor(lab_l=lab2[0], lab_a=lab2[1], lab_b=lab2[2])
    # delta_e = delta_e_cie1976(color1, color2)
    # print("de", delta_e)
    # assert abs(delta_e - ref) < 1.0e-12 * abs(delta_e)


def test_vector(tol=1.0e-14):
    rng = np.random.default_rng(0)
    lab1 = rng.random((3, 100))
    lab2 = rng.random((3, 100))
    refs = colorio.diff.cie76(lab1, lab2)
    for l1, l2, ref in zip(lab1.T, lab2.T, refs):
        val = colorio.diff.cie76(l1, l2)
        assert abs(val - ref) < tol * abs(ref)

    # test against reference
    norms = [65.33027889187598, 6.934455256389215, 1.1538228482108195]
    print(np.linalg.norm(refs, 1))
    print(np.linalg.norm(refs, 2))
    print(np.linalg.norm(refs, np.inf))
    assert abs(np.linalg.norm(refs, 1) - norms[0]) < tol * abs(norms[0])
    assert abs(np.linalg.norm(refs, 2) - norms[1]) < tol * abs(norms[1])
    assert abs(np.linalg.norm(refs, np.inf) - norms[2]) < tol * abs(norms[2])
