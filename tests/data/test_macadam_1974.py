import pytest

import colorio


def test_show():
    cs = colorio.cs.CIELAB
    plt = colorio.data.MacAdam1974().plot(cs)
    plt.show()
    plt.close()


@pytest.mark.parametrize(
    "cs_class,ref",
    [
        (colorio.cs.CAM16UCS, 18.862479231260807),
        (colorio.cs.CAM02UCS, 19.284192922142996),
        (colorio.cs.CAM02LCD, 19.126808698649114),
        (colorio.cs.OsaUcs, 19.847537952565744),
        (colorio.cs.CAM02SCD, 20.99363351524682),
        (colorio.cs.IPT, 21.269648149845857),
        (colorio.cs.JzAzBz, 23.211133853276525),
        (colorio.cs.RLAB, 23.524666960638214),
        (colorio.cs.SRLAB2, 23.85142864495729),
        (colorio.cs.CIELAB, 24.547740293433446),
        (colorio.cs.PROLAB, 27.02340205223324),
        (colorio.cs.CIELUV, 27.118673999331094),
        (colorio.cs.OKLAB, 32.71785472852726),
        (colorio.cs.ICtCp, 38.45053180241051),
        (colorio.cs.XYZ100, 50.86715667101449),
        (colorio.cs.XYY100, 84.30445664684683),
        (colorio.cs.CIEHCL, 86.73927353976319),
        (colorio.cs.CIELCH, 86.87548989979005),
    ],
)
def test_stress(cs_class, ref):
    res = colorio.data.MacAdam1974().stress(cs_class)
    print(cs_class)
    print(res)
    assert abs(res - ref) < 1.0e-13 * ref


if __name__ == "__main__":
    test_show()
    # test_residual()
