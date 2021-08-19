import json
import pathlib

import numpy as np

from ..helpers import ColorDistanceDataset


class Leeds(ColorDistanceDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "leeds.json") as f:
            data = json.load(f)

        xyz = np.asarray(data["xyz"])
        pairs = np.asarray(data["pairs"])

        super().__init__("Leeds", data["dv"], xyz[pairs])
