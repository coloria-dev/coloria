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

from ..helpers import ColorDistanceDataset


class BfdP(ColorDistanceDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "bfd-p.json") as f:
            data = json.load(f)

        # that's CIE 2-degree D64:
        self.whitepoint_xyz100 = data["reference_white"]

        xyz = np.asarray(data["xyz"])
        pairs = np.asarray(data["pairs"])

        super().__init__("BFD-P", data["dv"], xyz[pairs])
