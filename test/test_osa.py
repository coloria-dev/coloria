import numpy
import perfplot
import pytest

import colorio

numpy.random.seed(0)


@pytest.mark.parametrize(
    "xyz",
    [
        # difficult case that fails if the initial values aren't chosen carefully
        [12.0, 67.0, 20.0],
        100 * numpy.random.rand(3),
        100 * numpy.random.rand(3, 7),
        100 * numpy.random.rand(3, 4, 5),
    ],
)
def test_conversion(xyz):
    osa = colorio.OsaUcs()
    out = osa.to_xyz100(osa.from_xyz100(xyz))
    assert numpy.all(abs(xyz - out) < 1.0e-10)


def test_speed(N=2):
    numpy.random.seed(0)
    osa = colorio.OsaUcs()
    perfplot.plot(
        setup=lambda n: numpy.random.rand(3, n) * 100,
        equality_check=None,
        kernels=[
            osa.from_xyz100,
            numpy.cbrt,
        ],
        n_range=[2 ** n for n in range(N)],
        logx=True,
        # logy=True,
        relative_to=1
    )
    import tikzplotlib as tpl
    tpl.save("out.tex")


if __name__ == "__main__":
    test_speed(N=23)
