import json
import pathlib

import numpy as np

from ..helpers import HueLinearityDataset

this_dir = pathlib.Path(__file__).resolve().parent


class EbnerFairchild(HueLinearityDataset):
    def __init__(self):
        with open(this_dir / "ebner_fairchild.json") as f:
            data = json.load(f)

        arms = [
            np.column_stack([dat["reference xyz"], np.array(dat["same"]).T])
            for dat in data["data"]
        ]
        super().__init__("Ebner-Fairchild", data["white point"], arms)
