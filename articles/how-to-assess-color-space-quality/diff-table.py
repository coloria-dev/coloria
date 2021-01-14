# import matplotlib
import dufte
import matplotlib.pyplot as plt
import numpy
import tikzplotlib

import colorio

color_spaces = [
    colorio.cs.CAM02("UCS", 0.69, 20, 64 / numpy.pi / 5),
    colorio.cs.CAM16UCS(0.69, 20, 64 / numpy.pi / 5),
    colorio.cs.CIELAB(),
    colorio.cs.CIELUV(),
    colorio.cs.IPT(),
    colorio.cs.ICtCp(),
    colorio.cs.JzAzBz(),
    colorio.cs.OKLAB(),
    colorio.cs.OsaUcs(),
    colorio.cs.PROLAB(),
    colorio.cs.RLAB(),
    colorio.cs.XYY1(),
]

for cs in color_spaces:
    vals = [
        colorio.data.macadam_1942.stress(cs, 50),
        colorio.data.macadam_1974.stress(cs),
    ]
    print(f"{cs.name} & {vals[0]:.1f} & {vals[1]:.1f}\\\\")

labels = [cs.name for cs in color_spaces]
macadam_1942 = [colorio.data.macadam_1942.stress(cs, 50) for cs in color_spaces]
macadam_1974 = [colorio.data.macadam_1974.stress(cs) for cs in color_spaces]

plt.style.use(dufte.style)

x = numpy.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width / 2, macadam_1942, width, label="MacAdam 1942")
rects2 = ax.bar(x + width / 2, macadam_1974, width, label="MacAdam 1974")

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel("p_STRESS")
ax.set_xticks(x)
ax.set_xticklabels(labels)
plt.ylim(0, 100)
ax.legend()

fig.tight_layout()
# plt.show()
tikzplotlib.save("pstress.tex")
