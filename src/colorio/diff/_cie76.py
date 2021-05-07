import numpy as np


def cie76(lab1, lab2):
    diff = np.asarray(lab1) - np.asarray(lab2)
    return np.sqrt(np.einsum("i...,i...->...", diff, diff))
