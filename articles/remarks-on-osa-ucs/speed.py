import numpy as np
import perfplot

import colorio

np.random.seed(1)
osa = colorio.cs.OsaUcs()
cielab = colorio.cs.CIELAB()
# cam16 = colorio.CAM16(0.69, 20, L_A=64 / np.pi / 5)
ciecam02 = colorio.cs.CIECAM02(0.69, 20, L_A=64 / np.pi / 5)

perfplot.show(
    # Don't use np.random.rand(3, n) to avoid the CIECAM breakdown
    setup=lambda n: np.outer(np.random.rand(3) * 10, np.ones(n)),
    equality_check=None,
    kernels=[
        osa.to_xyz100,
        cielab.to_xyz100,
        # cam16.to_xyz100,
        lambda Jsh: ciecam02.to_xyz100(Jsh, "Jsh"),
        # np.cbrt,
    ],
    labels=["OSA-UCS", "CIELAB", "CIECAM02", "cbrt"],
    n_range=[2 ** n for n in range(23)],
    # relative_to=3
)
# import tikzplotlib as tpl
# tpl.save("out.tex")
