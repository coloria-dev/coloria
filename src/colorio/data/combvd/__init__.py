import json
import pathlib

import numpy as np

from ..helpers import ColorDistanceDataset


class COMBVD(ColorDistanceDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        components = [
            (this_dir / ".." / "bfd_p" / "bfd-p.json", 1.0),
            (this_dir / ".." / "leeds" / "leeds.json", 9.0),
            (this_dir / ".." / "rit_dupont" / "rit-dupont.json", 9.0),
            (this_dir / ".." / "witt" / "witt.json", 7.0),
        ]

        dist = []
        pairs = []
        weights = []
        for file_path, weight in components:
            with open(file_path) as f:
                data = json.load(f)
            dist.append(data["dv"])
            pairs.append(data["pairs"])
            weights.append(np.full(len(data["pairs"]), weight))

        dist = np.concatenate(dist)
        pairs = np.concatenate(pairs)
        weights = np.concatenate(weights)

        # CIECAM02 viewing conditions from the JzAzBz paper:
        self.L_A = 64
        self.c = 0.69
        self.Yb = 20

        super().__init__("COMBVD", dist, pairs, weights)
