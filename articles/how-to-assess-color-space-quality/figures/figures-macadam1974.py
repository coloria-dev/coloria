import matplotlib.pyplot as plt
import tikzplotlib

import colorio

data = [
    ("macadam1974-xyy.tex", colorio.cs.XYY(1)),
    ("macadam1974-cielab.tex", colorio.cs.CIELAB()),
    ("macadam1974-cam16.tex", colorio.cs.CAM16UCS(0.69, 20, 4.07)),
]
for filename, cs in data:
    colorio.data.MacAdam1974().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    tikzplotlib.save(
        filename,
        extra_axis_parameters=["width=0.25\\textwidth"],
    )
    plt.close()
