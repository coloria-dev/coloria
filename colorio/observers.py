import os

import numpy
import yaml


def cie_1931_2():
    """CIE 1931 standard observer, 2 degrees.
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, "data/observers/cie-1931-2.yaml")) as f:
        data = yaml.safe_load(f)
    data = numpy.array(data)
    return data[:, 0], data[:, 1:].T


def cie_1964_10():
    """CIE 1964 standard observer, 10 degrees.
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, "data/observers/cie-1964-10.yaml")) as f:
        data = yaml.safe_load(f)
    data = numpy.array(data)
    return data[:, 0], data[:, 1:].T
