"""
M.R. Luo, B. Rigg,
Chromaticity Discrimination Ellipses for Surface Colours,
Color Research and Application, Volume 11, Issue 1, Spring 1986, Pages 25-42,
<https://doi.org/10.1002/col.5080110107>.
"""
import json
import pathlib

import numpy as np

from ...cs import ColorCoordinates
from ..ellipse import EllipseDataset


class LuoRigg(EllipseDataset):
    def __init__(self, num_offset_points: int):
        # Extract ellipse centers and offsets from MacAdams data
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "luo-rigg.json") as f:
            data = json.load(f)

        xyy100_centers = []
        xyy100_points = []
        # collect the ellipse centers and offsets
        alpha = np.linspace(0.0, 2 * np.pi, num_offset_points, endpoint=False)
        circle_pts = np.array([np.cos(alpha), np.sin(alpha)])
        for data_set in data.values():
            # The set factor is the mean of the R values
            # set_factor = sum([dat[-1] for dat in data_set.values()]) / len(data_set)
            for x, y, Y, a, a_div_b, theta_deg, _ in data_set.values():
                theta = theta_deg * 2 * np.pi / 360
                a /= 1.0e4
                a *= (Y / 30) ** 0.2
                b = a / a_div_b
                # plot the ellipse
                c = np.array([x, y])
                J = np.array(
                    [
                        [+a * np.cos(theta), -b * np.sin(theta)],
                        [+a * np.sin(theta), +b * np.cos(theta)],
                    ]
                )
                offsets = np.dot(J, circle_pts)
                pts = (c + offsets.T).T

                xyy100_centers.append(np.array([*c, Y]))
                xyy100_points.append(np.array([*pts, np.full(pts.shape[1], Y)]).T)

        xyy100_centers = ColorCoordinates(np.array(xyy100_centers).T, "XYY100")
        xyy100_points = ColorCoordinates(np.array(xyy100_points).T, "XYY100")

        super().__init__("Luo-Rigg", xyy100_centers, xyy100_points)
