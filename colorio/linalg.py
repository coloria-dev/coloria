# -*- coding: utf-8 -*-
#
import numpy


def dot(A, x):
    '''Performs a dot product between an array of shape (m, n) and and array of
    shape (m, ...). The output has the same shape as the second argument.
    '''
    return numpy.einsum('ij,j...->i...', A, x)


def solve(A, x):
    '''Solves a linear equation system with a matrix of shape (n, n) and an
    array of shape (n, ...). The output has the same shape as the second
    argument.
    '''
    # https://stackoverflow.com/a/48387507/353337
    return numpy.linalg.solve(A, x.reshape(x.shape[0], -1)).reshape(x.shape)
