import matplotlib.pyplot as plt
import numpy as np

L_ = 25.0

t = np.linspace(4.1, 5.6, 100)
y = (L_ / 5.9 + 2 / 3 - t) ** 3 - 0.042**3 * (t**3 - 30)

plt.plot(t, y)
plt.gca().set_aspect("equal")
plt.grid()
plt.savefig("f.svg", transparent=True, bbox_inches="tight")
plt.show()
