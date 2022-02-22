"""
Mark D. Fairchild and Ping-Hsu Chen,
Brightness, Lightness, and Specifying Color in High-Dynamic-Range Scenes and Images
Munsell Color Science Laboratory, Chester F. Carlson Center for Imaging Science,
Rochester Institute of Technology, Rochester, NY, USA 14623-5604,
<https://doi.org/10.1117/12.872075>.
"""
import json
import pathlib
from typing import Type

import matplotlib.pyplot as plt
import numpy as np

from ...cs import ColorSpace
from ...illuminants import whitepoints_cie1931
from ..helpers import create_cs_class_instance, stress_absolute

this_dir = pathlib.Path(__file__).resolve().parent


class FairchildChen:
    def __init__(self, key: str):
        assert key in ["SL1", "SL2"]
        with open(this_dir / "fairchild_chen.json") as f:
            data = json.load(f)

        if key == "SL1":
            # From the article:
            #
            # > SL1 (4.8° viewing angle)
            #
            # > A projector was used to illuminate the configuration in a dark room and
            # > to modulate the luminances of the three patches in the center.
            #
            # CIECAM02 viewing conditions from the JzAzBz paper:
            self.c = 0.69  # "average"
            self.Y_b = 20
            self.L_A = 168
        else:
            assert key == "SL2"
            #
            # > SL2<100 & SL2>100 (2° viewing angle)
            #
            # > The luminance at the area of paper white was about 997 cd/m2.
            #
            # CIECAM02 viewing conditions from the JzAzBz paper:
            self.c = 0.69  # "average"
            self.Y_b = 20
            self.Lw = 997
            self.L_A = 199  # ~ Lw * Y_b / 100

        self.whitepoint_xyz100 = whitepoints_cie1931["D65"]

        self.data = data[key]
        self.key = key
        for key in self.data:
            self.data[key] = np.asarray(self.data[key])

    def plot(self, cs_class: Type[ColorSpace]):
        cs = create_cs_class_instance(
            cs_class, self.whitepoint_xyz100, self.c, self.Y_b, self.L_A
        )

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
        plt.legend()
        return plt

    def stress(self, cs_class: Type[ColorSpace]) -> float:
        cs = create_cs_class_instance(
            cs_class, self.whitepoint_xyz100, self.c, self.Y_b, self.L_A
        )

        # experimental lightness
        L = self.data["lightness"]
        # predicted lightness
        # Move L0 into origin for translation invariance
        assert cs.k0 is not None
        L0_ = cs.from_xyz100(np.zeros(3))[cs.k0]
        L_ = cs.from_xyz100(self.data["xyz"].T)[cs.k0] - L0_

        return stress_absolute(L, L_)
