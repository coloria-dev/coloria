import pathlib

import numpy as np
import yaml

this_dir = pathlib.Path(__file__).resolve().parent


def cie_1931_2():
    """CIE 1931 standard observer, 2 degrees."""
    with open(this_dir / "data/observers/cie-1931-2.yaml") as f:
        data = yaml.safe_load(f)
    data = np.array(data)
    return data[:, 0], data[:, 1:].T


def cie_1964_10():
    """CIE 1964 standard observer, 10 degrees."""
    with open(this_dir / "data/observers/cie-1964-10.yaml") as f:
        data = yaml.safe_load(f)
    data = np.array(data)
    return data[:, 0], data[:, 1:].T
