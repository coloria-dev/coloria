import numpy as np
import pytest
from scipy.optimize import dual_annealing, minimize

import colorio


def dot(a, b):
    """Take arrays `a` and `b` and form the dot product between the last axis
    of `a` and the first of `b`.
    """
    b = np.asarray(b)
    return np.dot(a, b.reshape(b.shape[0], -1)).reshape(a.shape[:-1] + b.shape[1:])


class TestLab(colorio.cs.ColorSpace):
    def __init__(self, x, whitepoint=colorio.illuminants.whitepoints_cie1931["D65"]):
        super().__init__("OPTLAB", ("L", "a", "b"), 0)
        self.whitepoint = whitepoint
        # self.p = x[0]
        # self.p = 1.0 / 3.0
        # self.M1 = x[0:9].reshape(3, 3)
        # self.M2 = x[9:18].reshape(3, 3)
        self.k0 = 0
        self.labels = ["L", "a", "b"]
        self.name = "TestLab"
        self.M = np.array([[0, 1, 0], [125 / 29, -125 / 29, 0], [0, 50 / 29, -50 / 29]])

    def from_xyz100(self, xyz):
        def f(t):
            delta = 6 / 29
            out = np.array(t, dtype=float)
            is_greater = out > delta**3
            out[is_greater] = 116 * np.cbrt(out[is_greater]) - 16
            out[~is_greater] = out[~is_greater] / (delta / 2) ** 3
            return out

        xyz = np.asarray(xyz)

        fxyz = f((xyz.T / self.whitepoint).T)
        return self.M @ fxyz

    # def to_xyz100(self, xyz):
    #     self.M1inv = np.linalg.inv(self.M1)
    #     self.M2inv = np.linalg.inv(self.M2)
    #     return dot(self.M1inv, dot(self.M2inv, xyz) ** (1.0 / self.p))


@pytest.mark.skip
def test_optimize(maxiter=1):
    # luo_rigg = colorio.data.LuoRigg(8)
    bfd_p = colorio.data.BfdP()
    leeds = colorio.data.Leeds()
    macadam_1942 = colorio.data.MacAdam1942(Y=50)
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

    def average(a, p):
        n = len(a)
        if p == 0:
            return np.prod(np.abs(a)) ** (1 / n)
        elif p == np.infty:
            return np.max(np.abs(a))
        return np.linalg.norm(a, p) / n ** (1 / p)

    def fun(x):
        cs = TestLab(x)
        # A typical hazard in this optimization is the collapse of the color space into
        # two dimensions (one lightness). This satisfies the hue linearity functionals
        # since all points are in one line then. The variant-a stress does not provide
        # enough counterweight. A brute-force remedy is to use small weights for the
        # hue-linearity functionals.
        variant = "a"
        d = np.array(
            [
                [1.0, bfd_p.stress(cs, variant)],
                [1.0, leeds.stress(cs, variant)],
                # [1.0, macadam_1942.stress(cs, variant)],
                [1.0, macadam_1974.stress(cs, variant)],
                [1.0, rit_dupont.stress(cs, variant)],
                [1.0, witt.stress(cs, variant)],
                #
                [0.1, average(ebner_fairchild.stress(cs), 2.0)],
                [0.1, average(hung_berns.stress(cs), 2.0)],
                [0.1, average(xiao.stress(cs), 2.0)],
                #
                [1.0, munsell.stress_lightness(cs)],
                [1.0, fairchild_chen.stress(cs)],
            ]
        )
        # make sure the weights add up to one
        d[:, 0] /= np.sum(d[:, 0])

        res = average(d[:, 0] * d[:, 1], 1.0)
        if np.isnan(res):
            res = 1.0e10
        # print(res)
        return res

    # x0 = np.array(
    #     [
    #         0.8189330101, 0.3618667424, -0.1288597137,
    #         0.0329845436, 0.9293118715, 0.0361456387,
    #         0.0482003018, 0.2643662691, 0.6338517070,
    #         #
    #         0.2104542553, +0.7936177850, -0.0040720468,
    #         +1.9779984951, -2.4285922050, +0.4505937099,
    #         +0.0259040371, +0.7827717662, -0.8086757660,
    #     ]
    # )

    # global search
    bounds = np.empty((18, 2))
    bounds[:, 0] = -3.0
    bounds[:, 1] = 3.0
    out = dual_annealing(fun, bounds)
    print("intermediate residual:")
    print(fun(out.x))
    # refine with bfgs
    out = minimize(
        fun,
        out.x,
        # method="Nelder-Mead",
        # method="Powell",
        # method="CG",
        method="BFGS",
        options={"maxiter": maxiter},
    )

    cs = TestLab(out.x)

    print()
    print("final residual:")
    print(fun(out.x))

    print()
    print("final values:")
    # print(out.x[0])
    print(out.x[0:9].reshape(3, 3))
    print(out.x[9:18].reshape(3, 3))

    print()
    print("final residuals:")
    d = {
        "  BFD-P............": bfd_p,
        "  Leeds............": leeds,
        "  MacAdam 1942.....": macadam_1942,
        "  MacAdam 1974.....": macadam_1974,
        "  RIT-DuPont.......": rit_dupont,
        "  Witt.............": witt,
    }
    for name, module in d.items():
        vals = [module.stress(cs, variant) for variant in ["a", "b", "c"]]
        string = name + " {:.3f}  {:.3f}  {:.3f}"
        print(string.format(*vals))
    print()
    d = {
        "Hung-Berns.......": hung_berns,
        "Ebner-Fairchild..": ebner_fairchild,
        "Xiao.............": xiao,
    }
    for name, module in d.items():
        stress = module.stress(cs)
        vals = [average(stress, p) for p in [1.0, 2.0, np.infty]]
        print("  " + name + " {:.3f}  {:.3f}  {:.3f}".format(*vals))
    print()
    print(f"  Munsell.......... {munsell.stress_lightness(cs):.3f}")
    print(f"  Fairchild-Chen... {fairchild_chen.stress(cs):.3f}")

    cs.show_primary_srgb_gradients()
    macadam_1942.show(cs)
    macadam_1974.show(cs)
    hung_berns.show(cs)
    ebner_fairchild.show(cs)
    xiao.show(cs)
    fairchild_chen.show(cs)


if __name__ == "__main__":
    test_optimize(maxiter=100000)
