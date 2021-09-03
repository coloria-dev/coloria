import json
import pathlib

import numpy as np

from ._helpers import SpectralData

this_dir = pathlib.Path(__file__).resolve().parent


# use 5nm as the default step size; it's consistent with the CIE standard
def cie_1931_2(stepsize: int = 5):
    return _from_file(this_dir / "data/observers/cie-1931-2.json", stepsize)


def cie_1964_10(stepsize: int = 5):
    return _from_file(this_dir / "data/observers/cie-1964-10.json", stepsize)


def _from_file(filename: pathlib.Path, stepsize: int):
    with open(filename) as f:
        data = json.load(f)

    lmbda_start, lmbda_end, lmbda_step = data["lambda_nm"]
    assert lmbda_step == 1
    lmbda = np.arange(lmbda_start, lmbda_end + 1, stepsize)
    vals = np.array(data["xyz"])[:, ::stepsize]

    return SpectralData(data["name"], lmbda, vals)
