import numpy as np
from numpy.typing import ArrayLike


class SpectralData:
    def __init__(self, name: str, lmbda_nm: ArrayLike, data: ArrayLike):
        self.name = name
        self.lmbda_nm = np.asarray(lmbda_nm)
        self.data = np.asarray(data)

    def __repr__(self) -> str:
        step = self.lmbda_nm[1] - self.lmbda_nm[0]
        return f"<{self.name} with {self.lmbda_nm[0]}:{step}:{self.lmbda_nm[-1]}nm>"
