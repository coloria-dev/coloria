import json
import pathlib

import numpy as np

from ._helpers import SpectralData

this_dir = pathlib.Path(__file__).resolve().parent


def cie_1931_2():
    with open(this_dir / "data/observers/cie-1931-2.json") as f:
        data = json.load(f)

    lmbda_start, lmbda_end, lmbda_step = data["lambda_nm"]
    lmbda = np.arange(lmbda_start, lmbda_end + 1, lmbda_step)
    return SpectralData(data["name"], lmbda, data["xyz"])


def cie_1964_10():
    with open(this_dir / "data/observers/cie-1964-10.json") as f:
        data = json.load(f)

    lmbda_start, lmbda_end, lmbda_step = data["lambda_nm"]
    lmbda = np.arange(lmbda_start, lmbda_end + 1, lmbda_step)
    return SpectralData(data["name"], lmbda, data["xyz"])
