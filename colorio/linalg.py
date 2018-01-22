# -*- coding: utf-8 -*-
#
import numpy


def dot(A, x):
    return numpy.einsum('ij,j...->i...', A, x)


def solve(A, x):
    return numpy.linalg.solve(A, x.reshape(x.shape[0], -1)).reshape(x.shape)
