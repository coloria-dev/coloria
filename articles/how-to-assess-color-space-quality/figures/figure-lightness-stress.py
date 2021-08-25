# import matplotlib
import dufte
import matplotlib.pyplot as plt
import numpy as np
import tikzplotlib

import colorio

ex = {
    "Fairchild-Chen (SL1) \\cite{fairchildchen}": colorio.data.FairchildChen("SL1"),
    "Fairchild-Chen (SL2)": colorio.data.FairchildChen("SL2"),
    "Munsell value": colorio.data.Munsell(),
}

cs = [
    colorio.cs.CAM02UCS,
    colorio.cs.CAM16UCS,
    colorio.cs.CIELAB,
    colorio.cs.CIELUV,
    # colorio.cs.ICtCp,
    colorio.cs.IPT,
    colorio.cs.JzAzBz,
    colorio.cs.OKLAB,
    colorio.cs.OsaUcs,
    colorio.cs.XYY1,
]


data_sets = {key: [data.stress(c) for c in cs] for key, data in ex.items()}


plt.style.use(dufte.style)

x = np.arange(len(cs))
n = len(data_sets)
bar_width = 0.8 / n

fig, ax = plt.subplots()
pos = np.linspace(-(n - 1) * bar_width / 2, (n - 1) * bar_width / 2, n, endpoint=True)
for (label, data), p in zip(data_sets.items(), pos):
    ax.bar(x + p, data, bar_width, label=label)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_title("l_{text{STRESS}}")
ax.set_xticks(x)
ax.set_xticklabels([c.name for c in cs])
plt.xticks(rotation=45, ha="right")
plt.xlim(-0.6, len(cs) - 1 + 0.6)
plt.ylim(0, 50)
ax.legend(framealpha=1, loc="upper left", bbox_to_anchor=(0, 1))

fig.tight_layout()
# plt.show()
tikzplotlib.save(
    "light_stress.tex",
    extra_axis_parameters=["width=\\textwidth", "height=0.5\\textwidth"],
    externalize_tables=True,
    override_externals=True,
    externals_search_path="./figures/",
)
