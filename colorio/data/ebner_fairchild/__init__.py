import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np

from ..helpers import Dataset, _compute_straight_line_stress, _plot_hue_linearity_data

this_dir = pathlib.Path(__file__).resolve().parent


class EbnerFairchild(Dataset):
    def __init__(self):
        with open(this_dir / "ebner_fairchild.json") as f:
            data = json.load(f)

        self.wp = np.array(data["white point"])
        self.d = [
            np.column_stack([dat["reference xyz"], np.array(dat["same"]).T])
            for dat in data["data"]
        ]

    def plot(self, cs):
        _plot_hue_linearity_data(self.d, self.wp, cs)
        plt.title(f"Ebner-Fairchild hue linearity data for {cs.name}")

    def stress(self, cs):
        return _compute_straight_line_stress(cs, self.wp, self.d)
