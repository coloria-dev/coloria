import matplotlib.pyplot as plt
import tikzplotlib

import colorio

data = [
    ("hung-berns-xyy.tex", colorio.cs.XYY1),
    ("hung-berns-cielab.tex", colorio.cs.CIELAB),
    ("hung-berns-oklab.tex", colorio.cs.OKLAB),
]
for filename, cs in data:
    colorio.data.HungBerns().plot(cs)
    plt.gca().set_aspect("equal")
    plt.title(cs.name)
    # plt.show()
    tikzplotlib.save(
        filename,
        extra_axis_parameters=["width=0.39\\textwidth", "ticks=none"],
        externalize_tables=True,
        override_externals=True,
        externals_search_path="./figures/",
    )
    plt.close()
