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

cs_labels = [
    "CAM02 (UCS)**",
    "CAM16 (UCS)**",
    "CIELAB*",
    "CIELUV*",
    # "$IC_tC_p$",
    "IPT",
    "$J_zA_zB_z$",
    "OKLAB",
    "OSA-UCS",
    "xyY",
]

data_sets = {
    key: [
        data.stress(
            colorio.cs.CAM02(
                "UCS",
                c=data.c,
                Y_b=data.Y_b,
                L_A=data.L_A,
                whitepoint=data.whitepoint_xyz100,
            )
        ),
        data.stress(
            colorio.cs.CAM16UCS(
                c=data.c, Y_b=data.Y_b, L_A=data.L_A, whitepoint=data.whitepoint_xyz100
            )
        ),
        data.stress(colorio.cs.CIELAB(whitepoint=data.whitepoint_xyz100)),
        data.stress(colorio.cs.CIELUV(whitepoint=data.whitepoint_xyz100)),
        # data.stress(colorio.cs.ICtCp()),
        data.stress(colorio.cs.IPT()),
        data.stress(colorio.cs.JzAzBz()),
        data.stress(colorio.cs.OKLAB()),
        data.stress(colorio.cs.OsaUcs()),
        data.stress(colorio.cs.XYY(1)),
    ]
    for key, data in ex.items()
}


plt.style.use(dufte.style)

x = np.arange(len(cs_labels))
n = len(data_sets)
bar_width = 0.8 / n

fig, ax = plt.subplots()
pos = np.linspace(-(n - 1) * bar_width / 2, (n - 1) * bar_width / 2, n, endpoint=True)
for (label, data), p in zip(data_sets.items(), pos):
    ax.bar(x + p, data, bar_width, label=label)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_title("l_{text{STRESS}}")
ax.set_xticks(x)
ax.set_xticklabels(cs_labels)
plt.xticks(rotation=45, ha="right")
plt.xlim(-0.6, len(cs_labels) - 1 + 0.6)
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
