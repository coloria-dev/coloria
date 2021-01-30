import numpy as np


def bisect(f, a, b, tol, max_num_steps=np.infty):
    fa = f(a)
    fb = f(b)

    assert np.all(np.logical_xor(fa > 0, fb > 0))
    is_fa_positive = fa > 0
    # sort points such that f(a) is negative, f(b) positive
    tmp = fa[is_fa_positive]
    fa[is_fa_positive] = fb[is_fa_positive]
    fb[is_fa_positive] = tmp
    #
    tmp = a[:, is_fa_positive]
    a[:, is_fa_positive] = b[:, is_fa_positive]
    b[:, is_fa_positive] = tmp

    diff = a - b
    dist2 = np.einsum("ij,ij->j", diff, diff)

    k = 0
    while True:
        if k >= max_num_steps:
            break

        c = (a + b) / 2
        fc = f(c)

        is_fc_positive = fc > 0
        a[:, ~is_fc_positive] = c[:, ~is_fc_positive]
        b[:, is_fc_positive] = c[:, is_fc_positive]

        diff = a - b
        dist2 = np.einsum("ij,ij->j", diff, diff)

        if np.all(dist2 < tol ** 2):
            break

        k += 1

    return a, b


def regula_falsi(f, a, b, tol, max_num_steps=np.infty):
    fa = f(a)
    fb = f(b)

    assert np.all(np.logical_xor(fa > 0, fb > 0))
    is_fa_positive = fa > 0
    # sort points such that f(a) is negative, f(b) positive
    tmp = fa[is_fa_positive]
    fa[is_fa_positive] = fb[is_fa_positive]
    fb[is_fa_positive] = tmp
    #
    tmp = a[:, is_fa_positive]
    a[:, is_fa_positive] = b[:, is_fa_positive]
    b[:, is_fa_positive] = tmp

    diff = a - b
    dist2 = np.einsum("ij,ij->j", diff, diff)

    k = 0
    while True:
        if k >= max_num_steps:
            break

        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)

        is_fc_positive = fc > 0
        a[:, ~is_fc_positive] = c[:, ~is_fc_positive]
        b[:, is_fc_positive] = c[:, is_fc_positive]

        diff = a - b
        dist2 = np.einsum("ij,ij->j", diff, diff)

        if np.all(dist2 < tol ** 2):
            break

        k += 1

    return a, b
