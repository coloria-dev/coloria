import numpy
import pytest

import colorio
import perfplot

numpy.random.seed(0)


@pytest.mark.parametrize(
    "xyz",
    [
        [0.0, 0.0, 0.0],
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
    assert numpy.all(abs(xyz - out) < 1.0e-9)


def test_speed(N=2):
    numpy.random.seed(1)
    osa = colorio.OsaUcs()
    cielab = colorio.CIELAB()
    # cam16 = colorio.CAM16(0.69, 20, L_A=64 / numpy.pi / 5)
    ciecam02 = colorio.CIECAM02(0.69, 20, L_A=64 / numpy.pi / 5)

    perfplot.plot(
        # Don't use numpy.random.rand(3, n) to avoid the CIECAM breakdown
        setup=lambda n: numpy.outer(numpy.random.rand(3), numpy.ones(n)) * 10,
        equality_check=None,
        kernels=[
            osa.to_xyz100,
            cielab.to_xyz100,
            # cam16.to_xyz100,
            lambda Jsh: ciecam02.to_xyz100(Jsh, "Jsh"),
            numpy.cbrt,
        ],
        labels=["OSA-UCS", "CIELAB", "CIECAM02", "cbrt"],
        n_range=[2 ** n for n in range(N)],
        logx=True,
        logy=True,
        # relative_to=3
    )
    # import tikzplotlib as tpl
    # tpl.save("out.tex")


if __name__ == "__main__":
    test_speed(N=2)
