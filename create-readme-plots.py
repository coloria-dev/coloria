import dufte
import matplotlib.pyplot as plt

import colorio

plt.style.use(dufte.style)


def normalize(string):
    return (
        string.replace("(", "")
        .replace(")", "")
        .replace(" ", "")
        .replace("_", "")
        .lower()
    )


data = [
    (colorio.cs.OKLAB(), 0.5, 2.0),
    (colorio.cs.CIELAB(), 50, 500),
    (colorio.cs.CAM16UCS(0.69, 20, 20), 50, 250),
]
for cs, lightness, camera_elevation in data:
    p = colorio.plot_rgb_slice(
        cs,
        lightness=lightness,
        camera_elevation=camera_elevation,
        n=101,
        off_screen=True,
    )
    # p.show()
    filename = f"srgb-gamut-slice-{normalize(cs.name)}.png"
    # tight? https://github.com/pyvista/pyvista-support/issues/487
    p.screenshot(filename, transparent_background=True)
    p.close()


data = [
    (colorio.cs.XYY(1), 0.5),
    (colorio.cs.OKLAB(), 0.5),
    (colorio.cs.JzAzBz(), 0.5),
]
for cs, lightness in data:
    plt = colorio.plot_visible_slice(cs, lightness=lightness)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    filename = f"visible-gamut-slice-{normalize(cs.name)}.png"
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    plt = colorio.data.Munsell().plot(cs, V=5)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    filename = f"munsell-{normalize(cs.name)}.svg"
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.MacAdam1942(50.0).plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    filename = f"macadam1942-{normalize(cs.name)}.svg"
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


for cs in [colorio.cs.XYY(1), colorio.cs.CIELAB(), colorio.cs.CAM16UCS(0.69, 20, 20)]:
    colorio.data.LuoRigg(8).plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    filename = f"luo-rigg-{normalize(cs.name)}.svg"
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.EbnerFairchild().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    filename = f"ebner-fairchild-{normalize(cs.name)}.svg"
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


for filename, cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.HungBerns().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    filename = f"hung-berns-{normalize(cs.name)}.svg"
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.Xiao().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    plt.savefig(f"xiao-{normalize(cs.name)}.svg", transparent=True, bbox_inches="tight")
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.MacAdam1974().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    filename = f"macadam1974-{normalize(cs.name)}.svg"
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


for filename, cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.FairchildChen("SL2").plot(cs)
    # plt.gca().set_aspect("equal")
    # plt.gca().grid(False)
    # plt.show()
    filename = f"fairchild-chen-{normalize(cs.name)}.svg"
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()
