import matplotlib.pyplot as plt
import tikzplotlib

import coloria

data = [
    ("macadam1974-xyy.tex", coloria.cs.XYY1),
    ("macadam1974-cielab.tex", coloria.cs.CIELAB),
    ("macadam1974-cam16.tex", coloria.cs.CAM16UCS),
]
for filename, cs in data:
    coloria.data.MacAdam1974().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(visible=False)
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
