from pathlib import Path

import matplotlib.pyplot as plt
import matplotx

import colorio
from colorio.cs import ColorCoordinates

plt.style.use(matplotx.styles.dufte)


def normalize(string):
    return (
        string.replace("(", "")
        .replace(")", "")
        .replace(" ", "")
        .replace("_", "")
        .lower()
    )


plot_dir = Path(__file__).resolve().parent / "plots"
plot_dir.mkdir(exist_ok=True)


data = [
    ("oklab", 0.5, 2.0),
    ("cielab", 50, 500),
    ("cam16ucs", 50, 250),
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
    # tight? https://github.com/pyvista/pyvista-support/issues/487
    p.screenshot(
        str(plot_dir / f"srgb-gamut-slice-{cs}.png"),
        transparent_background=True,
    )
    p.close()


data = [
    ("xyy100", 0.5),
    ("oklab", 0.5),
    ("jzazbz", 0.5),
]
for cs, lightness in data:
    plt = colorio.plot_visible_slice(cs, lightness=lightness)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    plt.savefig(
        plot_dir / f"visible-gamut-slice-{cs}.png",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    plt = colorio.data.Munsell().plot(cs, V=5)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    plt.savefig(
        plot_dir / f"munsell-{normalize(cs.name)}.svg",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.MacAdam1942(50.0).plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    plt.savefig(
        plot_dir / f"macadam1942-{normalize(cs.name)}.svg",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


for cs in ["xyy1", "cielab", "cam16ucs"]:
    colorio.data.LuoRigg(8).plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    plt.savefig(
        plot_dir / f"luo-rigg-{cs}.svg",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.EbnerFairchild().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    plt.savefig(
        plot_dir / f"ebner-fairchild-{normalize(cs.name)}.svg",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.HungBerns().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    plt.savefig(
        plot_dir / f"hung-berns-{normalize(cs.name)}.svg",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.Xiao().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    plt.savefig(
        plot_dir / f"xiao-{normalize(cs.name)}.svg",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.MacAdam1974().plot(cs)
    plt.gca().set_aspect("equal")
    plt.gca().grid(False)
    # plt.show()
    plt.savefig(
        plot_dir / f"macadam1974-{normalize(cs.name)}.svg",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


for cs in [colorio.cs.XYY1, colorio.cs.CIELAB, colorio.cs.CAM16UCS]:
    colorio.data.FairchildChen("SL2").plot(cs)
    # plt.gca().set_aspect("equal")
    # plt.gca().grid(False)
    # plt.show()
    plt.savefig(
        plot_dir / f"fairchild-chen-{normalize(cs.name)}.svg",
        transparent=True,
        bbox_inches="tight",
    )
    plt.close()


# illuminants
illu = colorio.illuminants.a()
plt.plot(illu.lmbda_nm, illu.data, label="A")
illu = colorio.illuminants.c()
plt.plot(illu.lmbda_nm, illu.data, label="C")
illu = colorio.illuminants.d50()
plt.plot(illu.lmbda_nm, illu.data, label="D50")
illu = colorio.illuminants.d65()
plt.plot(illu.lmbda_nm, illu.data, label="D65")
illu = colorio.illuminants.f2()
plt.plot(illu.lmbda_nm, illu.data, label="F2")
plt.xlabel("wavelength [nm]")
matplotx.line_labels()
plt.savefig(plot_dir / "illuminants.svg", transparent=True, bbox_inches="tight")
plt.close()


# observer
cols = [
    ColorCoordinates([30.0, 0.0, 0.0], "xyz100").convert("srgb-linear", mode="clip"),
    ColorCoordinates([0.0, 30.0, 0.0], "xyz100").convert("srgb-linear", mode="clip"),
    ColorCoordinates([0.0, 0.0, 30.0], "xyz100").convert("srgb-linear", mode="clip"),
]
obs = colorio.observers.cie_1931_2()
plt.plot(obs.lmbda_nm, obs.data[0], color=cols[0], label="$\\overline{x}$")
plt.plot(obs.lmbda_nm, obs.data[1], color=cols[1], label="$\\overline{y}$")
plt.plot(obs.lmbda_nm, obs.data[2], color=cols[2], label="$\\overline{z}$")
plt.xlabel("wavelength [nm]")
plt.legend()
plt.savefig(
    plot_dir / "cie-standard-observer-2.svg",
    transparent=True,
    bbox_inches="tight",
)
plt.show()
