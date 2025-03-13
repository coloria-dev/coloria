import matplotlib.pyplot as plt
import tikzplotlib

import coloria

data = [
    ("hung-berns-xyy.tex", coloria.cs.XYY1),
    ("hung-berns-cielab.tex", coloria.cs.CIELAB),
    ("hung-berns-oklab.tex", coloria.cs.OKLAB),
]
for filename, cs in data:
    coloria.data.HungBerns().plot(cs)
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
