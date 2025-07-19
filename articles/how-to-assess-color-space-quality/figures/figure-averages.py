import dufte
import matplotlib.pyplot as plt
import numpy as np
import tikzplotlib

plt.style.use(dufte.style)


def pnorm(a, p):
    n = len(a)
    if p == 0:
        return np.prod(np.abs(a)) ** (1 / n)
    return (np.sum(a**p) / n) ** (1 / p)


n = 20

p = np.linspace(-5.0, 5.0, 101)
# p = np.linspace(1.0, 10.0, 100)

a = np.zeros(n)
a[0] = 1.0
y1 = [pnorm(a, P) for P in p]

a = np.ones(n)
a[1] = 1.0e-4
y5 = [pnorm(a, P) for P in p]

a = np.linspace(0.0, 1.0, n, endpoint=True)
a[0] = 1.0e-4
y2 = [pnorm(a, P) for P in p]

a = np.full(n, 1.0)
y3 = [pnorm(a, P) for P in p]

a = np.full(n, 0.0)
y6 = [pnorm(a, P) for P in p]

plt.xlabel(r"\(p\)")
plt.plot(p, y3, label=r"\(s_{\text{equal1}}\)")
plt.plot(p, y5, label=r"\(s_{\text{outlier-small}}\)")
plt.plot(p, y2, label=r"\(s_{\text{uniform}}\)")
plt.plot(p, y1, label=r"\(s_{\text{outlier1}}\)")
plt.plot(p, y6, label=r"\(s_{\text{equal0}}\)")
plt.ylim(-0.05, 1.1)
# plt.gca().set_aspect("equal")

# xticks = list(plt.xticks()[0])
# print(xticks)
# plt.xticks(xticks + [1.0, 2.0], labels=["lol")
# plt.xticklabels(list(plt.xticklabels()[0]) + ["1 (avg)", "2 (RMS)"])


# dufte.legend()
plt.legend()
# plt.show()
tikzplotlib.save(
    "averages.tex",
    extra_axis_parameters=["width=\\textwidth", "height=0.4\\textwidth"],
    externalize_tables=True,
    override_externals=True,
    externals_search_path="./figures/",
)
