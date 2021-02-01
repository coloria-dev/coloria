import matplotlib.pyplot as plt
import dufte

import colorio

plt.style.use(dufte.style)

# cs = colorio.cs.CIELAB()
# colorio.save_rgb_slice("srgb-gamut-slice-cielab.png", cs, lightness=50, n=101)
#
# cs = colorio.cs.OKLAB()
# colorio.save_rgb_slice("srgb-gamut-slice-oklab.png", cs, lightness=0.5, n=101)
#
# cs = colorio.cs.CAM16UCS(0.69, 20, 4.07)
# colorio.save_rgb_slice("srgb-gamut-slice-cam16.png", cs, lightness=50, n=101)

# cs = colorio.cs.XYY(1)
# colorio.plot_rgb_slice(cs, lightness=0.5, n=101)
# colorio.plot_visible_slice(cs, lightness=0.5)
# plt.savefig("visible-gamut-slice-xyy.png", transparent=True, bbox_inches="tight")
#
# cs = colorio.cs.JzAzBz()
# colorio.plot_rgb_slice(cs, lightness=50, n=101)
# colorio.plot_visible_slice(cs, lightness=50)
# plt.savefig("visible-gamut-slice-jzazbz.png", transparent=True, bbox_inches="tight")

# cs = colorio.cs.OKLAB()
# colorio.plot_rgb_slice(cs, lightness=0.5, n=101)
# colorio.plot_visible_slice(cs, lightness=0.5)
# plt.savefig("visible-gamut-slice-oklab.png", transparent=True, bbox_inches="tight")

data = [
    ("munsell-xyy.png", colorio.cs.XYY(1)),
    ("munsell-cielab.png", colorio.cs.CIELAB()),
    ("munsell-cam16.png", colorio.cs.CAM16UCS(0.69, 20, 4.07)),
]
for filename, cs in data:
    colorio.data.Munsell().plot(cs, V=5)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()
