# import matplotlib
import dufte
import matplotlib.pyplot as plt
import numpy as np
import tikzplotlib

import colorio

bfdp = colorio.data.BfdP()
combvd = colorio.data.COMBVD()
leeds = colorio.data.Leeds()
macadam_1942 = colorio.data.MacAdam1942(Y=50.0)
macadam_1974 = colorio.data.MacAdam1974()
rit_dupont = colorio.data.RitDupont()
witt = colorio.data.Witt()

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
    colorio.cs.XYY(1),
]

xlabels = [cs.name for cs in color_spaces]
data_sets = {
    "BFD-P \\cite{luorigg}": [bfdp.stress(cs) for cs in color_spaces],
    "COMBVD": [combvd.stress(cs) for cs in color_spaces],
    "Leeds \\cite{leeds}": [leeds.stress(cs) for cs in color_spaces],
    "MacAdam \\cite{macadam1942} ($Y=50$)": [
        macadam_1942.stress(cs) for cs in color_spaces
    ],
    "MacAdam \\cite{macadam1974}": [macadam_1974.stress(cs) for cs in color_spaces],
    "RIT--DuPont \\cite{berns}": [rit_dupont.stress(cs) for cs in color_spaces],
    "Witt \\cite{witt}": [witt.stress(cs) for cs in color_spaces],
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
ax.set_title("$p_{STRESS}$")
ax.set_xticks(x)
ax.set_xticklabels(xlabels)
plt.xticks(rotation=45, ha="left")
plt.xlim(-0.6, len(xlabels) - 1 + 0.6)
plt.ylim(0, 100)
ax.legend()

fig.tight_layout()
# plt.show()
tikzplotlib.save(
    "pstress.tex",
    extra_axis_parameters=["width=\\textwidth", "height=0.6\\textwidth"],
    externalize_tables=True,
    override_externals=True,
    externals_search_path="./figures/",
)
