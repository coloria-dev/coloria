import matplotlib.pyplot as plt

import colorio

# cs = colorio.cs.CIELAB()
# colorio.save_rgb_slice("srgb-gamut-slice-cielab.png", cs, lightness=50, n=101)
#
# cs = colorio.cs.OKLAB()
# colorio.save_rgb_slice("srgb-gamut-slice-oklab.png", cs, lightness=0.5, n=101)
#
# cs = colorio.cs.CAM16UCS(0.69, 20, 4.07)
# colorio.save_rgb_slice("srgb-gamut-slice-cam16.png", cs, lightness=50, n=101)

cs = colorio.cs.CIELAB()
colorio.plot_rgb_slice(cs, lightness=50, n=101)
colorio.plot_visible_slice(cs, lightness=50)
plt.savefig("visible-gamut-slice-cielab.png", transparent=True, bbox_inches="tight")