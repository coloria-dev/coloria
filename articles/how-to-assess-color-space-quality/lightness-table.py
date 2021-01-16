# import matplotlib
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
    colorio.cs.IPT(),
    colorio.cs.ICtCp(),
    colorio.cs.JzAzBz(),
    colorio.cs.OKLAB(),
    colorio.cs.OsaUcs(),
    colorio.cs.PROLAB(),
    colorio.cs.RLAB(),
    colorio.cs.XYY(1),
]

# for cs in color_spaces:
#     vals = [
#         colorio.data.fairchild_chen.stress(cs, "SL1"),
#         colorio.data.fairchild_chen.stress(cs, "SL2"),
#     ]
#     print(f"{cs.name} & {vals[0]:.1f} & {vals[1]:.1f}\\\\")

xlabels = [cs.name for cs in color_spaces]
data_sets = {
    "Fairchild-Chen (SL1) \\cite{fairchildchen}": [
        colorio.data.fairchild_chen.stress(cs, "SL1") for cs in color_spaces
    ],
    "Fairchild-Chen (SL2)": [
        colorio.data.fairchild_chen.stress(cs, "SL2") for cs in color_spaces
    ],
    "Munsell value": [colorio.data.munsell.stress_lightness(cs) for cs in color_spaces],
}

plt.style.use(dufte.style)

x = np.arange(len(xlabels))
n = len(data_sets)
bar_width = 0.8 / n

fig, ax = plt.subplots()
pos = np.linspace(-(n - 1) * bar_width / 2, (n - 1) * bar_width / 2, n, endpoint=True)
for (label, data), p in zip(data_sets.items(), pos):
    ax.bar(x + p, data, bar_width, label=label)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel("l_STRESS")
ax.set_xticks(x)
ax.set_xticklabels(xlabels)
plt.xticks(rotation=45)
plt.ylim(0)
ax.legend()

fig.tight_layout()
# plt.show()
tikzplotlib.save(
    "light_stress.tex",
    extra_axis_parameters=["width=\\textwidth", "height=0.5\\textwidth"],
)
