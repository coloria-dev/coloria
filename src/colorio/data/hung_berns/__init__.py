import json
import pathlib

import numpy as np

from ...illuminants import whitepoints_cie1931
from ..helpers import HueLinearityDataset


class HungBerns(HueLinearityDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent
        with open(this_dir / "table3.json") as f:
            data = json.load(f)

        # CIECAM02 viewing conditions from the JzAzBz paper:
        self.L_A = 10
        self.c = 0.525
        self.Yb = 20

        arms = [np.array(list(color.values())).T for color in data.values()]
        super().__init__("Hung-Berns", whitepoints_cie1931["C"], arms)
