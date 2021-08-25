import matplotlib.pyplot as plt
import tikzplotlib

import colorio

data = [
    ("macadam1974-xyy.tex", colorio.cs.XYY1),
    ("macadam1974-cielab.tex", colorio.cs.CIELAB),
    ("macadam1974-cam16.tex", colorio.cs.CAM16UCS),
]
for filename, cs in data:
    colorio.data.MacAdam1974().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    plt.title(cs.name)
    tikzplotlib.save(
        filename,
        extra_axis_parameters=["width=0.4\\textwidth", "ticks=none"],
        extra_lines_start=["\\scriptsize"],
        externalize_tables=True,
        override_externals=True,
        externals_search_path="./figures/",
    )
    plt.close()
