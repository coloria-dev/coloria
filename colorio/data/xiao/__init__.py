import json
import pathlib

import numpy as np

from ..helpers import HueLinearityDataset


class Xiao(HueLinearityDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "averages.json") as f:
            data = json.load(f)

        neutral_gray = np.array(data.pop("neutral-gray")[0])
        arms = np.moveaxis(list(data.values()), 1, 2)
        super().__init__("Xiao", neutral_gray, arms)
