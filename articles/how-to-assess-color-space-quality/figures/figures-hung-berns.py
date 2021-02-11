import matplotlib.pyplot as plt
import tikzplotlib

import colorio

data = [
    ("hung-berns-xyy.tex", colorio.cs.XYY(1)),
    ("hung-berns-cielab.tex", colorio.cs.CIELAB()),
    ("hung-berns-oklab.tex", colorio.cs.OKLAB()),
]
for filename, cs in data:
    colorio.data.HungBerns().plot(cs)
    plt.gca().set_aspect("equal")
    # plt.show()
    tikzplotlib.save(
        filename,
        extra_axis_parameters=["width=0.33\\textwidth"],
    )
    plt.close()
