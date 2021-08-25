import colorio


def test_stress():
    data = colorio.data.COMBVD()
    cs = colorio.cs.CIELAB
    ref = 43.92861966265739
    res = data.stress(cs)
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


def test_stress_relative():
    data = colorio.data.COMBVD()
    cs = colorio.cs.CIELAB

    ref = 44.370141147103546
    res = data.stress(cs, variant="relative")
    print(res)
    assert abs(res - ref) < 1.0e-14 * ref


if __name__ == "__main__":
    # test_show()
    test_stress()
