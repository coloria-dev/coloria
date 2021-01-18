import pathlib

import matplotlib.pyplot as plt
import numpy as np
import yaml

from ..helpers import Dataset, _compute_straight_line_stress, _plot_hue_linearity_data


class Xiao(Dataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        data = []
        with open(this_dir / "averages.yaml") as f:
            data = yaml.safe_load(f)
        data = np.array(list(data.values()))

        # Use Xiao's 'neutral gray' as white point.
        with open(this_dir / "neutral_gray.yaml") as f:
            ng_data = np.array(yaml.safe_load(f))
        self.ng = np.sum(ng_data, axis=0) / np.prod(ng_data.shape[:1])
        self.data = np.moveaxis(data, 1, 2)

    def plot(self, cs):
        _plot_hue_linearity_data(self.data, self.ng, cs)
        plt.title(f"Xiao hue linearity data for {cs.name}")

    def stress(self, cs):
        return _compute_straight_line_stress(cs, self.ng, self.data)
