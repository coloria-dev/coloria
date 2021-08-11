import dufte
import matplotlib.pyplot as plt
import numpy as np
import tikzplotlib

import colorio

diff_formulas = {
    "CIE76": colorio.diff.cie76,
    "CIE94": colorio.diff.cie94,
    "CIEDE2000": colorio.diff.ciede2000,
    "CMC 2:1": colorio.diff.cmc,
}

bfdp = colorio.data.BfdP()
combvd = colorio.data.COMBVD()
leeds = colorio.data.Leeds()
macadam_1942 = colorio.data.MacAdam1942(Y=50.0)
macadam_1974 = colorio.data.MacAdam1974()
rit_dupont = colorio.data.RitDupont()
witt = colorio.data.Witt()

xlabels = list(diff_formulas.keys())
data_sets = {
    "BFD-P \\cite{luorigg}": [
        bfdp.stress_lab_diff(cs) for cs in diff_formulas.values()
    ],
    "COMBVD": [combvd.stress_lab_diff(cs) for cs in diff_formulas.values()],
    "Leeds \\cite{leeds}": [leeds.stress_lab_diff(cs) for cs in diff_formulas.values()],
    "MacAdam \\cite{macadam1942} ($Y=50$)": [
        macadam_1942.stress_lab_diff(cs) for cs in diff_formulas.values()
    ],
    "MacAdam \\cite{macadam1974}": [
        macadam_1974.stress_lab_diff(cs) for cs in diff_formulas.values()
    ],
    "RIT--DuPont \\cite{berns}": [
        rit_dupont.stress_lab_diff(cs) for cs in diff_formulas.values()
    ],
    "Witt \\cite{witt}": [witt.stress_lab_diff(cs) for cs in diff_formulas.values()],
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
    "pstress-diff.tex",
    extra_axis_parameters=["width=\\textwidth", "height=0.6\\textwidth"],
    externalize_tables=True,
    override_externals=True,
    externals_search_path="./figures/",
)
