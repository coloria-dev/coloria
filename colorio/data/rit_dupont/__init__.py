import json
import pathlib

from ..helpers import ColorDifferenceDataset


class RitDupont(ColorDifferenceDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "rit-dupont.json") as f:
            data = json.load(f)

        super().__init__("RIT-DuPont", data["dv"], data["pairs"])
