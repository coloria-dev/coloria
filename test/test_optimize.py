import numpy as np
from scipy.optimize import minimize

import colorio


def dot(a, b):
    """Take arrays `a` and `b` and form the dot product between the last axis
    of `a` and the first of `b`.
    """
    b = np.asarray(b)
    return np.dot(a, b.reshape(b.shape[0], -1)).reshape(a.shape[:-1] + b.shape[1:])


class TestLab:
    def __init__(self, x):
        self.p = x[0]
        self.M1 = x[1:10].reshape(3, 3)
        self.M2 = x[10:19].reshape(3, 3)
        self.k0 = 0
        self.labels = ["L", "a", "b"]
        self.name = "TestLab"

    def from_xyz100(self, xyz):
        return dot(self.M2, dot(self.M1, xyz) ** self.p)


def test_optimize(maxiter=1):
    luo_rigg = colorio.data.LuoRigg(8)
    macadam_1942 = colorio.data.MacAdam1942()
    macadam_1974 = colorio.data.MacAdam1974()
    rit_dupont = colorio.data.RitDupont()
    witt = colorio.data.Witt()
    #
    ebner_fairchild = colorio.data.EbnerFairchild()
    hung_berns = colorio.data.HungBerns()
    xiao = colorio.data.Xiao()
    #
    munsell = colorio.data.Munsell()
    fairchild_chen = colorio.data.FairchildChen("SL2")

    def fun(x):
        # res = np.average(colorio.data.ebner_fairchild.stress(TestLab()))
        # res = colorio.data.macadam_1942.stress(TestLab(), 50)
        cs = TestLab(x)
        res = (
            # luo_rigg.stress(cs)
            macadam_1942.stress(cs, 50)
            # + macadam_1974.stress(cs)
            # + witt.stress(cs)
            # + rit_dupont.stress(cs)
            #
            # + np.average(ebner_fairchild.stress(cs))
            # + np.average(hung_berns.stress(cs))
            # + np.average(xiao.stress(cs))
            #
            # + munsell.stress_lightness(cs)
            # + fairchild_chen.stress(cs)
        )
        if np.isnan(res):
            res = 1.0e10
        # print()
        print(res)
        # print(cs.p)
        # print(cs.M1)
        # print(cs.M2)
        return res

    # np.random.seed(1)
    x0 = np.random.rand(19)

    # x0 = np.array(
    #     [
    #         1.0 / 3.0,
    #         #
    #         0.8189330101,
    #         0.3618667424,
    #         -0.1288597137,
    #         0.0329845436,
    #         0.9293118715,
    #         0.0361456387,
    #         0.0482003018,
    #         0.2643662691,
    #         0.6338517070,
    #         #
    #         0.2104542553,
    #         +0.7936177850,
    #         -0.0040720468,
    #         +1.9779984951,
    #         -2.4285922050,
    #         +0.4505937099,
    #         +0.0259040371,
    #         +0.7827717662,
    #         -0.8086757660,
    #     ]
    # )
    # x0 = np.array(
    #     [
    #         1.0,
    #         #
    #         1.0, 0.0, 0.0,
    #         0.0, 1.0, 0.0,
    #         0.0, 0.0, 1.0,
    #         #
    #         1.0, 0.0, 0.0,
    #         0.0, 1.0, 0.0,
    #         0.0, 0.0, 1.0,
    #     ]
    # )

    # print(fun(x0))

    out = minimize(
        fun,
        x0,
        # method="Nelder-Mead",
        # method="Powell",
        # method="CG",
        method="BFGS",
        options={"maxiter": maxiter},
    )
    print(out)

    cs = TestLab(out.x)

    print()
    print("final residual:")
    print(fun(out.x))

    print()
    print("final values:")
    print(out.x[0])
    print(out.x[1:10].reshape(3, 3))
    print(out.x[10:19].reshape(3, 3))

    print()
    print("final residuals:")
    print(luo_rigg.stress(cs))
    print(macadam_1942.stress(cs, 50))
    print(macadam_1974.stress(cs))
    print(rit_dupont.stress(cs))
    print(witt.stress(cs))
    print()
    print(np.average(hung_berns.stress(cs)))
    print(np.average(ebner_fairchild.stress(cs)))
    print(np.average(xiao.stress(cs)))
    print()
    print(munsell.stress_lightness(cs))
    print(fairchild_chen.stress(cs))

    luo_rigg.show(cs)
    macadam_1942.show(cs)
    macadam_1974.show(cs)
    hung_berns.show(cs)
    ebner_fairchild.show(cs)
    xiao.show(cs)
    fairchild_chen.show(cs)


if __name__ == "__main__":
    test_optimize(maxiter=10000)
