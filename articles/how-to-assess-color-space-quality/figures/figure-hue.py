import dufte
import matplotlib.pyplot as plt
import numpy as np
import tikzplotlib

import colorio

ex = {
    "Hung--Berns \\cite{hung}": colorio.data.HungBerns(),
    "Ebner--Fairchild \\cite{ebner}": colorio.data.EbnerFairchild(),
    "Xiao et al. \\cite{xiao}": colorio.data.Xiao(),
}

cs_labels = [
    "CAM02 (UCS)**",
    "CAM16 (UCS)**",
    "CIELAB*",
    "CIELUV*",
    "$IC_tC_p$",
    "IPT",
    "$J_zA_zB_z$",
    "OKLAB",
    "OSA-UCS",
    "xyY",
]

data_sets = {
    key: [
        data.stress(colorio.cs.CAM02UCS),
        data.stress(colorio.cs.CAM16UCS),
        data.stress(colorio.cs.CIELAB),
        data.stress(colorio.cs.CIELUV),
        data.stress(colorio.cs.ICtCp),
        data.stress(colorio.cs.IPT),
        data.stress(colorio.cs.JzAzBz),
        data.stress(colorio.cs.OKLAB),
        data.stress(colorio.cs.OsaUcs),
        data.stress(colorio.cs.XYY1),
    ]
    for key, data in ex.items()
}


plt.style.use(dufte.style)

# the label locations:
x = np.arange(len(cs_labels))
n = len(data_sets)
bar_width = 0.8 / n

fig, ax = plt.subplots()

pos = np.linspace(-(n - 1) * bar_width / 2, (n - 1) * bar_width / 2, n, endpoint=True)
# cat 20 color pairs
color_pairs = [
    ("#1f77b4", "#aec7e8"),
    ("#ff7f0e", "#ffbb78"),
    ("#2ca02c", "#98df8a"),
    ("#d62728", "#ff9896"),
    ("#9467bd", "#c5b0d5"),
    ("#8c564b", "#c49c94"),
    ("#e377c2", "#f7b6d2"),
    ("#7f7f7f", "#c7c7c7"),
    ("#bcbd22", "#dbdb8d"),
    ("#17becf", "#9edae5"),
]

for (label, data), p, cols in zip(data_sets.items(), pos, color_pairs):
    average = np.array([np.average(item) for item in data])
    maxval = np.array([np.max(item) for item in data])
    ax.bar(x + p, average, bar_width, label=label, color=cols[0], zorder=5)
    ax.bar(
        x + p,
        maxval - average,
        bar_width,
        bottom=average,
        # label=label,
        color=cols[1],
        zorder=5,
    )

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_title("$h_{text{STRESS}}$")
ax.yaxis.set_label_coords(-0.1, 1.02)
plt.xticks(x, rotation=45, ha="right")
ax.set_xticklabels(cs_labels)
plt.xlim(-0.6, len(cs_labels) - 1 + 0.6)
plt.ylim(0, 20)
ax.legend(framealpha=1, loc="upper right", bbox_to_anchor=(1, 1))

plt.gcf().tight_layout()
# plt.show()
tikzplotlib.save(
    "hstress.tex",
    extra_axis_parameters=["width=\\textwidth", "height=0.5\\textwidth"],
    externalize_tables=True,
    override_externals=True,
    externals_search_path="./figures/",
)
