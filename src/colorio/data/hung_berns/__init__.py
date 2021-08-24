"""
Po-Chieh Hung, Roy S. Berns,
Determination of Constant Hue Loci for a CRT Gamut and Their Predictions Using
Color Appearance Spaces,
Color Research and Application, Volume 20, Issue 5 October 1995, Pages 285-295,
<https://doi.org/10.1002/col.5080200506>.
"""
import json
import pathlib

import numpy as np

from ...illuminants import whitepoints_cie1931
from ..hue_linearity import HueLinearityDataset


class HungBerns(HueLinearityDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent
        with open(this_dir / "table3.json") as f:
            data = json.load(f)

        # From the article:
        #
        # > A BARCO CalibratorTM 19" high-resolution color monitor driven by a SUN-3/260
        # > workstation and a SUN TAAC-I graphics accelerator was used as the stimulus
        # > display.
        #
        # > All observations were made in a darkened room.
        #
        # Those are the CIECAM02 viewing conditions as given in the JzAzBz paper:
        self.L_A = 10
        self.c = 0.525  # "dark"
        self.Y_b = 20

        arms = [np.array(list(color.values())).T for color in data.values()]
        super().__init__("Hung-Berns", whitepoints_cie1931["C"], arms)
