import json
import pathlib

import numpy as np

from ..color_distance import ColorDistanceDataset


class Leeds(ColorDistanceDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "leeds.json") as f:
            data = json.load(f)

        xyz = np.asarray(data["xyz"])
        pairs = np.asarray(data["pairs"])

        self.whitepoint_xyz100 = data["reference_white"]

        # parameters as used in
        #
        # Melgosa, Huertas, Berns,
        # Performance of recent advanced color-difference formulas using the
        # standardized residual sum of squares index
        # J. Opt. Soc. Am. A/ Vol. 25, No. 7/ July 2008,
        # <https://doi.org/10.1364/JOSAA.25.001828>.
        #
        self.c = 0.69  # average
        self.L_A = 100
        self.Y_b = 20

        super().__init__("Leeds", data["dv"], xyz[pairs])
