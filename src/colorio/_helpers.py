import numpy as np


class SpectralData:
    def __init__(self, name, lmbda_start, lmbda_end, lmbda_step, xyz):
        self.name = name
        self.lmbda_nm = np.arange(lmbda_start, lmbda_end + 1, lmbda_step)
        self.xyz = np.asarray(xyz)

    def __repr__(self):
        lmbda_nm_step = self.lmbda_nm[1] - self.lmbda_nm[0]
        return f"<{self.name} with {lmbda_nm_step}nm resolution>"
