import dufte
import matplotlib.pyplot as plt

x = [0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000]
y1 = [
    0.000386412008083425,
    0.0575982900045346,
    0.124436895988765,
    0.191784891998395,
    0.255525363012566,
    0.320106088009197,
    0.381167707993882,
    0.447597951002535,
    0.50900621600158,
    0.585287615002017,
    0.649983358001919,
]
y2 = [
    0.000280210006167181,
    0.0540978200006066,
    0.116504425997846,
    0.1806468910072,
    0.242063231999055,
    0.300529030006146,
    0.365901857003337,
    0.430936441989616,
    0.485103309008991,
    0.555676405012491,
    0.616663155989954,
]

plt.style.use(dufte.style)

plt.plot(x, y1, label="without fixes")
plt.plot(x, y2, label="with fixes")
plt.xlabel("Number of input samples")
dufte.ylabel("Runtime [s]")
dufte.legend()

plt.savefig("cam16-fixes-speed-comparison.svg", transparent=True, bbox_inches="tight")
plt.show()
