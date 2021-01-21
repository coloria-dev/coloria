import numpy as np


def dot(a, b):
    """Take arrays `a` and `b` and form the dot product between the last axis
    of `a` and the first of `b`.
    """
    b = np.asarray(b)
    return np.dot(a, b.reshape(b.shape[0], -1)).reshape(a.shape[:-1] + b.shape[1:])


def solve(A, x):
    """Solves a linear equation system with a matrix of shape (n, n) and an
    array of shape (n, ...). The output has the same shape as the second
    argument.
    """
    # https://stackoverflow.com/a/48387507/353337
    x = np.asarray(x)
    return np.linalg.solve(A, x.reshape(x.shape[0], -1)).reshape(x.shape)


def divide_zero(a, b):
    """Division, but set to 0 if b == 0.
    """
    return np.divide(a, b, out=np.zeros_like(a), where=b != 0.0)
