import json
import pathlib

import numpy as np

from ..helpers import HueLinearityDataset


class EbnerFairchild(HueLinearityDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent
        with open(this_dir / "ebner_fairchild.json") as f:
            data = json.load(f)

        arms = [
            np.column_stack([dat["reference xyz"], np.array(dat["same"]).T])
            for dat in data["data"]
        ]

        # CIECAM02 viewing conditions from the JzAzBz paper:
        self.L_A = 14
        self.c = 0.525
        self.Yb = 20

        super().__init__("Ebner-Fairchild", data["white point"], arms)
