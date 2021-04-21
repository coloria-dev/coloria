import numpy as np
import pytest

import colorio


@pytest.mark.parametrize(
    "lab1, lab2, ref",
    [
        ((50.0000, 2.6772, -79.7751), (50.0000, 0.0000, -82.7485), 1.738736105726153),
        ((50.0000, -1.1848, -84.8006), (50.0000, 0.0000, -82.7485), 0.8833292683548156),
        ((50.0000, -1.0000, 2.0000), (50.0000, 0.0000, 0.0000), 2.8793003178961936),
        ((50.0000, 2.4900, -0.0010), (50.0000, -2.4900, 0.0011), 6.578380802119273),
        ((50.0000, 2.5000, 0.0000), (56.0000, -27.0000, -3.0000), 38.06176941087411),
        ((50.0000, 2.5000, 0.0000), (58.0000, 24.0000, 15.0000), 33.334208357392015),
        ((50.0000, 2.5000, 0.0000), (50.0000, 3.2592, 0.3350), 1.0533532771968366),
        ((61.2901, 3.7196, -5.3901), (61.4292, 2.2480, -4.9620), 1.7655589559175304),
        ((90.9257, -0.5406, -0.9208), (88.6381, -0.8985, -0.7239), 0.9900703671254522),
        ((6.7747, -0.2908, -2.4247), (5.8714, -0.0985, -2.2286), 0.9528008883107042),
        ((2.0776, 0.0795, -1.1350), (0.9033, -0.0636, -0.5514), 1.427772909322429),
    ],
)
def test(lab1, lab2, ref):
    print(lab1, lab2)
    print(ref)
    val = colorio.diff.cmc(lab1, lab2)
    print(val)
    assert abs(val - ref) < 1.0e-13 * abs(ref)

    # from colormath.color_diff import delta_e_cmc
    # from colormath.color_objects import LabColor

    # color1 = LabColor(*lab1)
    # color2 = LabColor(*lab2)
    # delta_e = delta_e_cmc(color1, color2)
    # print("de", delta_e)
    # assert abs(delta_e - ref) < 1.0e-12 * abs(delta_e)


def test_vector():
    rng = np.random.default_rng(0)
    lab1 = rng.random((3, 100))
    lab2 = rng.random((3, 100))
    refs = colorio.diff.cmc(lab1, lab2)
    for l1, l2, ref in zip(lab1.T, lab2.T, refs):
        val = colorio.diff.cmc(l1, l2)
        assert abs(val - ref) < 1.0e-14 * abs(ref)

    # test against reference
    norms = [85.42009136248436, 9.083007601505594, 1.5307485776155307]
    print(np.linalg.norm(refs, 1))
    print(np.linalg.norm(refs, 2))
    print(np.linalg.norm(refs, np.inf))
    assert abs(np.linalg.norm(refs, 1) - norms[0]) < 1.0e-14 * abs(norms[0])
    assert abs(np.linalg.norm(refs, 2) - norms[1]) < 1.0e-14 * abs(norms[1])
    assert abs(np.linalg.norm(refs, np.inf) - norms[2]) < 1.0e-14 * abs(norms[2])
