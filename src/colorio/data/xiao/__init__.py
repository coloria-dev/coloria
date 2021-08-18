"""
Kaida Xiao, Sophie Wuerger, Chenyang Fu, Dimosthenis Karatzas,
Unique Hue Data for Colour Appearance Models. Part I: Loci of Unique Hues and Hue
Uniformity,
Color Reseach and Application, Volume 36, Issue 5, October 2011, Pages 316-323,
<https://doi.org/10.1002/col.20637>
"""
import json
import pathlib

import numpy as np

from ..helpers import HueLinearityDataset


class Xiao(HueLinearityDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "averages.json") as f:
            data = json.load(f)

        whitepoint_xyz100 = np.array(data["neutral-gray"][0])
        whitepoint_xyz100 /= whitepoint_xyz100[1]
        whitepoint_xyz100 *= 100

        # The whitepoint as given in the article, table 3.
        # Note that this does not equal the "neutral gray" as given in the data,
        # [57.74156756756757, 59.86000000000005, 77.54005405405404],
        # scaled [96.46102166, 100.0 , 129.53567333] and
        # also slighty diverges from the value given in [Safdar et al., Perceptually
        # uniform color space for image signals including high dynamic range and wide
        # gamut] (the JzAzBz paper). The value [97.313, 100, 138.596] is used there.
        whitepoint_xyz100 = np.array([98.0, 100.0, 139.7])
        # CIECAM02 surround parameters, as given in the article
        self.Lw = 114.6
        self.Yb = 20
        self.surrounding = 0.59  # "dim"
        self.L_A = 23  # ~ Lw * Yb / 100

        data.pop("neutral-gray", None)
        arms = np.moveaxis(list(data.values()), 1, 2)
        super().__init__("Xiao", whitepoint_xyz100, arms)
