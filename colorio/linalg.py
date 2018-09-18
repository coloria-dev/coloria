# -*- coding: utf-8 -*-
#
import numpy


def dot(a, b):
    """Take arrays `a` and `b` and form the dot product between the last axis
    of `a` and the first of `b`.
    """
    b = numpy.asarray(b)
    return numpy.dot(a, b.reshape(b.shape[0], -1)).reshape(a.shape[:-1] + b.shape[1:])


def solve(A, x):
    """Solves a linear equation system with a matrix of shape (n, n) and an
    array of shape (n, ...). The output has the same shape as the second
    argument.
    """
    # https://stackoverflow.com/a/48387507/353337
    x = numpy.asarray(x)
    return numpy.linalg.solve(A, x.reshape(x.shape[0], -1)).reshape(x.shape)
