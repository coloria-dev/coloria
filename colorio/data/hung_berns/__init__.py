import pathlib

import numpy as np
import yaml

from ...illuminants import whitepoints_cie1931
from ..helpers import Dataset, _compute_straight_line_stress, _plot_hue_linearity_data


class HungBerns(Dataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent
        with open(this_dir / "table3.yaml") as f:
            data = yaml.safe_load(f)

        self.wp = whitepoints_cie1931["C"]
        self.d = [np.array(list(color.values())).T for color in data.values()]

    def plot(self, cs):
        _plot_hue_linearity_data(self.d, self.wp, cs)

    def stress(self, cs):
        return _compute_straight_line_stress(cs, self.wp, self.d)
