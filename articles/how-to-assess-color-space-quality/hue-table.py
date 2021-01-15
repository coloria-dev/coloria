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
    colorio.cs.RLAB(),
    colorio.cs.XYY(1),
]

for cs in color_spaces:
    vals = [
        colorio.data.hung_berns.stress(cs),
        colorio.data.ebner_fairchild.stress(cs),
        colorio.data.xiao.stress(cs),
    ]
    avg = np.average(np.concatenate(vals))
    vals = [np.average(val) for val in vals]
    # print(f"{cs.name}    {sum(vals)}")
    print(f"{cs.name} & {vals[0]:.1f} & {vals[1]:.1f} & {vals[2]:.1f} & {avg:.1f}\\\\")

labels = [cs.name for cs in color_spaces]
data_sets = {
    "Hung-Berns": [
        np.average(colorio.data.hung_berns.stress(cs)) for cs in color_spaces
    ],
    "Ebner-Fairchild": [
        np.average(colorio.data.ebner_fairchild.stress(cs)) for cs in color_spaces
    ],
    "Xiao et al.": [np.average(colorio.data.xiao.stress(cs)) for cs in color_spaces],
}

plt.style.use(dufte.style)

# the label locations:
x = np.arange(len(labels))
n = len(data_sets)
bar_width = 0.8 / n

fig, ax = plt.subplots()

pos = np.linspace(-(n - 1) * bar_width / 2, (n - 1) * bar_width / 2, n, endpoint=True)
for (label, data), p in zip(data_sets.items(), pos):
    ax.bar(x + p, data, bar_width, label=label)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel("l_STRESS")
ax.set_xticks(x)
ax.set_xticklabels(labels)
# iplt.ylim(0, 100)
ax.legend()

plt.gcf().tight_layout()
# plt.show()
tikzplotlib.save("lstress.tex")
