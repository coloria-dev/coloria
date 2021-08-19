import json
import pathlib

import numpy as np

from ..helpers import ColorDistanceDataset


class RitDupont(ColorDistanceDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "rit-dupont.json") as f:
            data = json.load(f)

        xyz = np.asarray(data["xyz"])
        pairs = np.asarray(data["pairs"])

        super().__init__("RIT-DuPont", data["dv"], xyz[pairs])
