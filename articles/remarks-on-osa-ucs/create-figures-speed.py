import matplotlib.pyplot as plt
import numpy as np
import perfplot

import colorio

# import tikzplotlib


rng = np.random.default_rng(1)
osa = colorio.cs.OsaUcs()
cielab = colorio.cs.CIELAB()
# cam16 = colorio.CAM16(0.69, 20, L_A=64 / np.pi / 5)
ciecam02 = colorio.cs.CIECAM02(0.69, 20, L_A=64 / np.pi / 5)

b = perfplot.bench(
    # Don't use np.random.rand(3, n) to avoid the CIECAM breakdown
    setup=lambda n: np.outer(rng.random(3) * 10, np.ones(n)),
    equality_check=None,
    kernels=[
        osa.to_xyz100,
        cielab.to_xyz100,
        # cam16.to_xyz100,
        lambda Jsh: ciecam02.to_xyz100(Jsh, "Jsh"),
        np.cbrt,
    ],
    labels=["OSA-UCS", "CIELAB", "CIECAM02", "cbrt"],
    n_range=[2**n for n in range(23)],
)
b.plot()
plt.title("Runtime [s]")
plt.ylabel("")
plt.savefig("speed-absolute.svg", transparent=True, bbox_inches="tight")
# tikzplotlib.save(
#     "figures/speed-absolute.tex",
#     externalize_tables=True,
#     override_externals=True,
#     externals_search_path="./figures/",
# )
plt.close()


b.plot(relative_to=3)
plt.title("Runtime relative to cbrt")
plt.ylabel("")
plt.savefig("speed-relative.svg", transparent=True, bbox_inches="tight")
# tikzplotlib.save(
#     "figures/speed-relative.tex",
#     externalize_tables=True,
#     override_externals=True,
#     externals_search_path="./figures/",
# )
plt.close()
