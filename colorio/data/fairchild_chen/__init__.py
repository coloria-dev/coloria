import pathlib

import dufte
import matplotlib.pyplot as plt
import numpy as np
import yaml


def load():
    this_dir = pathlib.Path(__file__).resolve().parent
    with open(this_dir / "fairchild_chen.yaml") as f:
        data = yaml.safe_load(f)

    for key in data:
        data[key]["lightness"] = np.asarray(data[key]["lightness"])
        data[key]["xyz"] = np.asarray(data[key]["xyz"])
    return data


def show(*args):
    plt.figure()
    plot(*args)
    plt.show()
    plt.close()


def savefig(filename, *args):
    plt.figure()
    plot(*args)
    plt.savefig(filename, transparent=True, bbox_inches="tight")
    plt.close()


def plot(cs, key: str):
    assert key in ["SL1", "SL2"]
    plt.style.use(dufte.style)

    data = load()[key]

    # experimental lightness
    L = data["lightness"]
    # predicted lightness
    L_ = cs.from_xyz100(data["xyz"].T)[cs.k0]

    Y = data["xyz"][:, 1]
    plt.plot(Y, L, "o", label="experimental")

    alpha = np.dot(L, L_) / np.dot(L, L)
    plt.plot(Y, L_ / alpha, "-", label=f"{cs.name} (scaled)")

    plt.xlabel("Y")
    plt.ylabel("lightness")
    plt.title(f"Fairchild-Chen {key} lightness data")
    # dufte.legend()
    plt.legend()
    plt.show()


def residual(cs, key: str):
    d = load()[key]

    # experimental lightness
    L = d["lightness"]
    # predicted lightness
    # Move L0 into origin for translation invariance
    L0_ = cs.from_xyz100(np.zeros(3))[cs.k0]
    L_ = cs.from_xyz100(d["xyz"].T)[cs.k0] - L0_

    alpha = np.dot(L, L_) / np.dot(L, L)
    diff = alpha * L - L_
    return np.sqrt(np.dot(diff, diff) / np.dot(L_, L_))


def stress(*args):
    return 100 * residual(*args)
