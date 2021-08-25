import dufte
import matplotlib.pyplot as plt
import tikzplotlib

import colorio

plt.style.use(dufte.style)


for cs in [colorio.cs.IPT, colorio.cs.CIELAB, colorio.cs.OsaUcs]:
    colorio.data.Munsell().plot_lightness(cs)
    plt.title(cs.name)
    tikzplotlib.clean_figure()
    tikzplotlib.save(
        f"munsell-lightness-{cs.name.lower()}.tex",
        extra_axis_parameters=["width=0.33\\textwidth", "height=0.3\\textwidth"],
        extra_lines_start=["\\scriptsize"],
        externalize_tables=True,
        override_externals=True,
        externals_search_path="./figures/",
    )
    plt.close()
