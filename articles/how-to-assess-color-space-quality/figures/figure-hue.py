import dufte
import matplotlib.pyplot as plt
import numpy as np
import tikzplotlib

import colorio

color_spaces = [
    colorio.cs.CAM02("UCS", 0.69, 20, 64 / np.pi / 5),
    colorio.cs.CAM16UCS(0.69, 20, 64 / np.pi / 5),
    colorio.cs.CIELAB(),
    colorio.cs.CIELUV(),
    colorio.cs.ICtCp(),
    colorio.cs.IPT(),
    colorio.cs.JzAzBz(),
    colorio.cs.OKLAB(),
    colorio.cs.OsaUcs(),
    colorio.cs.PROLAB(),
    colorio.cs.XYY(1),
]

# for cs in color_spaces:
#     vals = [
#         colorio.data.hung_berns.stress(cs),
#         colorio.data.ebner_fairchild.stress(cs),
#         colorio.data.xiao.stress(cs),
#     ]
#     avg = np.average(np.concatenate(vals))
#     vals = [np.average(val) for val in vals]
#     # print(f"{cs.name}    {sum(vals)}")
#     print(f"{cs.name} & {vals[0]:.1f} & {vals[1]:.1f} & {vals[2]:.1f} & {avg:.1f}\\\\")

hung_berns = colorio.data.HungBerns()
ebner_fairchild = colorio.data.EbnerFairchild()
xiao = colorio.data.Xiao()

labels = [cs.name for cs in color_spaces]
data_sets = {
    "Hung--Berns \\cite{hung}": [hung_berns.stress(cs) for cs in color_spaces],
    "Ebner--Fairchild \\cite{ebner}": [
        ebner_fairchild.stress(cs) for cs in color_spaces
    ],
    "Xiao et al. \\cite{xiao}": [xiao.stress(cs) for cs in color_spaces],
}

plt.style.use(dufte.style)

# the label locations:
x = np.arange(len(labels))
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
    minval = np.array([np.min(item) for item in data])
    ax.bar(
        x + p,
        average,
        bar_width,
        label=label,
        color=cols[0],
        zorder=5
    )
    ax.bar(
        x + p,
        maxval - average,
        bar_width,
        bottom=average,
        # label=label,
        color=cols[1],
        zorder=5
    )

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_title("$h_{STRESS}$")
ax.yaxis.set_label_coords(-0.1, 1.02)
plt.xticks(x, rotation=45, ha="right")
ax.set_xticklabels(labels)
plt.xlim(-0.6, len(labels) - 1 + 0.6)
plt.ylim(0, 25)
ax.legend()

plt.gcf().tight_layout()
# plt.show()
tikzplotlib.save(
    "hstress.tex", extra_axis_parameters=["width=\\textwidth", "height=0.5\\textwidth"]
)
