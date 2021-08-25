"""
F. Ebner and M.D. Fairchild,
Finding Constant Hue Surfaces in Color Space,
SPIE/IS&T Electronic Imaging, (1998),
<https://doi.org/10.1117/12.298269>.
"""
import json
import pathlib

import numpy as np

from ..hue_linearity import HueLinearityDataset


class EbnerFairchild(HueLinearityDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent
        with open(this_dir / "ebner_fairchild.json") as f:
            data = json.load(f)

        arms = [
            np.column_stack([dat["reference xyz"], np.array(dat["same"]).T])
            for dat in data["data"]
        ]

        # From the article:
        #
        # > The experiment was performed on a Sun Sparc Ultra 1 with a 20" GDM-20E20
        # > monitor. The resolution was 1152x900 with a refresh rate of 76 Hz and a
        # > 0.31mm-phosphor trio pitch. The CRT display was kept on during the extent of
        # > the experiment (which ran from May through July) to minimize drift from
        # > power cycling. The layout of the interface is shown in figure 1. The
        # > stimulus squares for both test and reference subtended an angle of 4.2
        # > degrees.
        # > The white border had a luminance of 71 cd/m2. The luminance factor of the
        # > background gray was .35 (25 cd/m2).
        #
        # > The CRT was set up to have a white point near D65, with a luminance of
        # > 71 cd/m2. The experiment took place under dark surround conditions. The
        # > walls of the room were covered with black felt to eliminate possibility of
        # > reflection from the CRT.
        #
        self.L_A = 25  # ~ 71 * Y_b / 100
        self.c = 0.525  # "dark"
        self.Y_b = 35

        super().__init__("Ebner-Fairchild", data["white point"], arms)
