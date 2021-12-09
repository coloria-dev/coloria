from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


class SpectralData:
    def __init__(self, lmbda_nm: ArrayLike, data: ArrayLike, name: str | None = None):
        self.name = "spectral data" if name is None else name
        self.lmbda_nm = np.asarray(lmbda_nm)
        self.data = np.asarray(data)

        n = len(self.lmbda_nm)
        assert n == self.data.shape[-1], (
            f"Data mismatch for {name}: "
            f"len(lmbda) = {n}, data.shape[-1] == {self.data.shape[-1]}"
        )

    def __repr__(self) -> str:
        step = self.lmbda_nm[1] - self.lmbda_nm[0]
        return f"<{self.name} with {self.lmbda_nm[0]}:{step}:{self.lmbda_nm[-1]}nm>"
