import numpy


def cie76(lab1, lab2):
    diff = numpy.asarray(lab1) - numpy.asarray(lab2)
    return numpy.sqrt(numpy.einsum("i...,i...->...", diff, diff))
