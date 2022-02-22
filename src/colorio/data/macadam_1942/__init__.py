"""
David L. MacAdam,
Visual Sensitivities to Color Differences in Daylight,
Journal of the Optical Society of America,
Volume 32, May, 1942, Number 5,
https://doi.org/10.1364/JOSA.32.000247
"""
import json
import pathlib
from typing import Type

import matplotlib.pyplot as plt
import numpy as np

from ...cs import ColorCoordinates, ColorSpace, convert
from ...illuminants import whitepoints_cie1931
from ..color_distance import ColorDistanceDataset
from ..ellipse import _plot_ellipses
from ..helpers import create_cs_class_instance


class MacAdam1942(ColorDistanceDataset):
    def __init__(self, Y):
        self.Y = Y

        # Whitepoint info from
        #
        # Muhammad Safdar, Guihua Cui, Youn Jin Kim, and Ming Ronnier Luo,
        # Perceptually uniform color space for image signals including high dynamic
        # range and wide gamut,
        # Optics Express Vol. 25, Issue 13, pp. 15131-15151 (2017),
        # <https://doi.org/10.1364/OE.25.015131>.
        #
        self.whitepoint_xyz100 = whitepoints_cie1931["C"]

        # CIECAM02 viewing conditions from the JzAzBz paper:
        self.L_A = 24
        self.c = 0.69
        self.Y_b = 20

        # Extract ellipse centers and offsets from MacAdams data
        this_dir = pathlib.Path(__file__).resolve().parent
        with open(this_dir / "table3.json") as f:
            data = json.load(f)
        #
        # collect the ellipse centers and offsets
        self.xy_centers = []
        self.xy_offsets = []
        for datak in data:
            # collect ellipse points
            _, _, _, _, delta_y_delta_x, delta_s = np.array(datak["data"]).T
            offset = (
                np.array([np.ones(delta_y_delta_x.shape[0]), delta_y_delta_x])
                / np.sqrt(1 + delta_y_delta_x**2)
                * delta_s
            )
            if offset.shape[1] < 2:
                continue
            self.xy_centers.append(np.array([datak["x"], datak["y"]]))
            self.xy_offsets.append(np.column_stack([+offset, -offset]))

        # create pairs
        xyy_pairs = []
        for center, offset in zip(self.xy_centers, self.xy_offsets):
            for pt in center + offset.T:
                xyy_pairs.append([[*center, Y], [*pt, Y]])
        xyy_pairs = ColorCoordinates(np.array(xyy_pairs).T, "xyy100")

        xyz_pairs = convert(xyy_pairs, "xyz100")

        data = xyz_pairs.data.T
        super().__init__("MacAdam (1942)", np.ones(len(data)), data)

    def plot(self, cs_class: Type[ColorSpace], ellipse_scaling: float = 10.0):
        cs = create_cs_class_instance(
            cs_class, self.whitepoint_xyz100, self.c, self.Y_b, self.L_A
        )

        Y = self.Y
        centers_points = []
        for c, off in zip(self.xy_centers, self.xy_offsets):
            ctr = np.array([*c, Y])
            p = (c + off.T).T
            pts = np.array([*p, np.full(p.shape[1], Y)])
            centers_points.append(
                ColorCoordinates(np.column_stack([ctr, pts]), "xyy100")
            )

        _plot_ellipses(cs, centers_points, ellipse_scaling)
        plt.title(f"MacAdam ellipses for {cs.name}")

        # cs.plot_visible_slice(
        #     lightness,
        #     outline_prec=outline_prec,
        #     fill_color=visible_gamut_fill_color,
        # )
        # if plot_srgb_gamut:
        #     cs.plot_rgb_slice(lightness)
        return plt
