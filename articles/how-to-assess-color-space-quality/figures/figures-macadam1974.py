import matplotlib.pyplot as plt
import tikzplotlib

import colorio

data = [
    ("macadam1974-xyy.tex", "XYY", colorio.cs.XYY(1)),
    ("macadam1974-cielab.tex", "CIELAB", colorio.cs.CIELAB()),
    ("macadam1974-cam16.tex", "CAM16", colorio.cs.CAM16UCS(0.69, 20, 4.07)),
]
for filename, title, cs in data:
    colorio.data.MacAdam1974().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    plt.title(title)
    # turn off ticks to save space

    # plt.show()
    tikzplotlib.save(
        filename,
        extra_axis_parameters=["width=0.4\\textwidth", "ticks=none"],
        extra_lines_start=["\\scriptsize"],
        externalize_tables=True,
        override_externals=True,
        externals_search_path="./figures/",
    )
    plt.close()

# filename, cs = "macadam1974-xyy.tex", colorio.cs.XYY(1)
# colorio.data.MacAdam1974().plot(cs)
# plt.gca().set_aspect("equal")
# plt.gca().grid(False)
# plt.show()
# plt.show()
# tikzplotlib.save(
#     filename,
#     extra_axis_parameters=["width=0.25\\textwidth"],
#     extra_lines_start=["\\scriptsize"],
#     externalize_tables=True,
#     override_externals=True,
#     externals_search_path="./figures/",
# )
# plt.close()
