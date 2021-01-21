import json
import pathlib

from ..helpers import ColorDifferenceDataset


class Leeds(ColorDifferenceDataset):
    def __init__(self):
        this_dir = pathlib.Path(__file__).resolve().parent

        with open(this_dir / "leeds.json") as f:
            data = json.load(f)

        super().__init__("Leeds", data["dv"], data["pairs"])
