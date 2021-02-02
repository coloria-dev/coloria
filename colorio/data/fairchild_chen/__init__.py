import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np

from ..helpers import Dataset

this_dir = pathlib.Path(__file__).resolve().parent


class FairchildChen(Dataset):
    def __init__(self, key: str):
        assert key in ["SL1", "SL2"]
        with open(this_dir / "fairchild_chen.json") as f:
            data = json.load(f)

        self.data = data[key]
        self.key = key
        for key in self.data:
            self.data[key] = np.asarray(self.data[key])

    def plot(self, cs):
        # experimental lightness
        L = self.data["lightness"]
        # predicted lightness
        L_ = cs.from_xyz100(self.data["xyz"].T)[cs.k0]

        alpha = np.dot(L, L_) / np.dot(L, L)
        Y = self.data["xyz"][:, 1]
        plt.plot(Y, L * alpha, "o", label="experimental (scaled)")

        plt.plot(Y, L_, "-", label=f"{cs.name}")

        plt.xlabel("Y")
        plt.ylabel(cs.labels[cs.k0], rotation=0)
        plt.title(f"Fairchild-Chen {self.key} lightness data")
        # dufte.legend()
        plt.legend()

    def stress(self, cs):
        # experimental lightness
        L = self.data["lightness"]
        # predicted lightness
        # Move L0 into origin for translation invariance
        L0_ = cs.from_xyz100(np.zeros(3))[cs.k0]
        L_ = cs.from_xyz100(self.data["xyz"].T)[cs.k0] - L0_

        alpha = np.dot(L, L_) / np.dot(L, L)
        diff = alpha * L - L_
        return 100 * np.sqrt(np.dot(diff, diff) / np.dot(L_, L_))
