import json
import pathlib

from ..helpers import ColorDifferenceDataset


class BfdP(ColorDifferenceDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "bfd-p.json") as f:
            data = json.load(f)

        super().__init__("BFD-P", data["dv"], data["pairs"])
