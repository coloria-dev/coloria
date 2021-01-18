import pathlib

import matplotlib.pyplot as plt
import numpy as np
import yaml

from ..helpers import Dataset, _compute_straight_line_stress, _plot_hue_linearity_data


class Xiao(Dataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "averages.yaml") as f:
            data = yaml.safe_load(f)

        self.ng = np.array(data.pop("neutral-gray")[0])
        data = np.array(list(data.values()))
        self.data = np.moveaxis(data, 1, 2)

    def plot(self, cs):
        _plot_hue_linearity_data(self.data, self.ng, cs)
        plt.title(f"Xiao hue linearity data for {cs.name}")

    def stress(self, cs):
        return _compute_straight_line_stress(cs, self.ng, self.data)
