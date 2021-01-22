import numpy as np
import dufte
import matplotlib.pyplot as plt
import tikzplotlib

plt.style.use(dufte.style)


def pnorm(a, p):
    n = len(a)
    if p == 0:
        return np.prod(np.abs(a)) ** (1 / n)
    return (np.sum(a ** p) / n) ** (1 / p)


np.random.seed(0)

n = 20

p = np.linspace(0.0, 5.0, 100)
# p = np.linspace(1.0, 10.0, 100)

a = np.zeros(n)
a[0] = 1.0
y1 = [pnorm(a, P) for P in p]

a = np.linspace(0.0, 1.0, n, endpoint=True)
y2 = [pnorm(a, P) for P in p]

a = np.full(n, 1.0)
y3 = [pnorm(a, P) for P in p]

plt.xlabel("$p$")
plt.plot(p, y3, label="$s_{equal}$")
plt.plot(p, y2, label="$s_{uniform}$")
plt.plot(p, y1, label="$s_{outlier}$")
plt.ylim(0.0, 1.1)
# plt.gca().set_aspect("equal")

# xticks = list(plt.xticks()[0])
# print(xticks)
# plt.xticks(xticks + [1.0, 2.0], labels=["lol")
# plt.xticklabels(list(plt.xticklabels()[0]) + ["1 (avg)", "2 (RMS)"])


# dufte.legend()
plt.legend()
# plt.show()
tikzplotlib.save(
   "norms.tex", extra_axis_parameters=["width=\\textwidth", "height=0.4\\textwidth"]
   )
