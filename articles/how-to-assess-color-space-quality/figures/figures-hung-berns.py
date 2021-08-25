import matplotlib.pyplot as plt
import tikzplotlib

import colorio

data = [
    ("hung-berns-xyy.tex", "xyY", colorio.cs.XYY1),
    ("hung-berns-cielab.tex", "CIELAB", colorio.cs.CIELAB),
    ("hung-berns-oklab.tex", "OKLAB", colorio.cs.OKLAB),
]
for filename, title, cs in data:
    colorio.data.HungBerns().plot(cs)
    plt.gca().set_aspect("equal")
    plt.title(title)
    # plt.show()
    tikzplotlib.save(
        filename,
        extra_axis_parameters=["width=0.39\\textwidth", "ticks=none"],
        externalize_tables=True,
        override_externals=True,
        externals_search_path="./figures/",
    )
    plt.close()
