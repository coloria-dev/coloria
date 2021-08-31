import json
import pathlib

from ._helpers import SpectralData

this_dir = pathlib.Path(__file__).resolve().parent


def cie_1931_2():
    with open(this_dir / "data/observers/cie-1931-2.json") as f:
        data = json.load(f)
    return SpectralData(data["name"], *data["lambda_nm"], data["xyz"])


def cie_1964_10():
    with open(this_dir / "data/observers/cie-1964-10.json") as f:
        data = json.load(f)
    return SpectralData(data["name"], *data["lambda_nm"], data["xyz"])
