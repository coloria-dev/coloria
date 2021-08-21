"""
C. Alder, K.P. Chaing, T.F. Chong, E. Coates, A.A. Khalili and B. Rigg,
Uniform Chromaticity Scales - New Experimental Data,
Journal of The Society of Dyers and Colourists,
Volume 98 January 1982,
<https://doi.org/10.1111/J.1478-4408.1982.TB03602.X>
"""
import json
import pathlib

import numpy as np

from ..color_distance import ColorDistanceDataset


class BfdP(ColorDistanceDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "bfd-p.json") as f:
            data = json.load(f)

        # surround parameters as used in
        #
        # Melgosa, Huertas, Berns,
        # Performance of recent advanced color-difference formulas using the
        # standardized residual sum of squares index
        # J. Opt. Soc. Am. A/ Vol. 25, No. 7/ July 2008,
        # <https://doi.org/10.1364/JOSAA.25.001828>.
        #
        self.c = 0.69  # average
        self.L_A = 100
        self.Yb = 20

        xyz = np.asarray(data["xyz"])
        pairs = np.asarray(data["pairs"])

        super().__init__("BFD-P", data["dv"], xyz[pairs])
