import dufte
import matplotlib.pyplot as plt
import numpy as np
import tikzplotlib

plt.style.use(dufte.style)


ymax = 2.0

x0 = np.linspace(0.0, 5.0, 101)
y0 = (1 - x0) ** 2

idx = y0 < ymax * 1.5
x0 = x0[idx]
y0 = y0[idx]

x1 = np.linspace(0.0, 5.0, 101)[1:]
y1 = (1 - x1) ** 2 / x1

idx = y1 < ymax * 1.5
x1 = x1[idx]
y1 = y1[idx]

x2 = np.linspace(0.0, 5.0, 101)
y2 = (1 - x2) ** 2 / x2**2

idx = y2 < ymax * 1.5
x2 = x2[idx]
y2 = y2[idx]

plt.plot(x0, y0, label="$(1-x)^2$")
plt.plot(x1, y1, label="$(1-x)^2 / x$")
plt.plot(x2, y2, label="$(1-x)^2 / x^2$")
plt.legend()
plt.ylim(0.0, ymax)
plt.gca().set_aspect("equal")
# plt.show()

tikzplotlib.save(
    "norm-scaling.tex",
    extra_axis_parameters=["width=\\textwidth", "height=0.4\\textwidth"],
    externalize_tables=True,
    override_externals=True,
    externals_search_path="./figures/",
)
