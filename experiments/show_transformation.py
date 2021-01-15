import argparse

import matplotlib.pyplot as plt
import meshzoo
import numpy as np
from dolfin import Function, FunctionSpace, Mesh, MeshEditor

import colorio


def _main():
    args = _parse_cmd_arguments()

    content = np.load(args.infile)

    data = content.item()["data"]
    n = content.item()["n"]

    # # plot statistics
    # axes0 = problem.get_ellipse_axes(alpha0).T.flatten()
    # plt.plot(axes0, label='axes lengths before')
    # axes1 = problem.get_ellipse_axes(out.x).T.flatten()
    # plt.plot(axes1, label='axes lengths opt')
    # plt.legend()
    # plt.grid()

    # Plot unperturbed MacAdam
    # colorio.plot_luo_rigg(
    #     ellipse_scaling=1,
    colorio.save_macadam(
        "macadam-native.png", ellipse_scaling=10, plot_rgb_triangle=False, n=n
    )

    points, cells = meshzoo.triangle(
        corners=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]), n=n
    )

    # https://bitbucket.org/fenics-project/dolfin/issues/845/initialize-mesh-from-vertices
    editor = MeshEditor()
    mesh = Mesh()
    editor.open(mesh, "triangle", 2, 2)
    editor.init_vertices(points.shape[0])
    editor.init_cells(cells.shape[0])
    for k, point in enumerate(points):
        editor.add_vertex(k, point[:2])
    for k, cell in enumerate(cells):
        editor.add_cell(k, cell)
    editor.close()

    V = FunctionSpace(mesh, "CG", 1)

    def get_u(alpha):
        n = V.dim()
        ax = alpha[:n]
        ay = alpha[n:]

        ux = Function(V)
        ux.vector().set_local(ax)
        ux.vector().apply("")

        uy = Function(V)
        uy.vector().set_local(ay)
        uy.vector().apply("")
        return ux, uy

    # Plot perturbed MacAdam
    def transform(XY, data=data):
        is_solo = len(XY.shape) == 1
        if is_solo:
            XY = np.array([XY]).T
        # print(XY)
        ux, uy = get_u(data)
        out = np.array([[ux(x, y) for x, y in XY.T], [uy(x, y) for x, y in XY.T]])
        if is_solo:
            out = out[..., 0]
        return out

    # colorio.plot_luo_rigg(
    #     ellipse_scaling=1,
    plt.figure()
    colorio.plot_macadam(
        ellipse_scaling=10,
        # xy_to_2d=problem.pade2d.eval,
        xy_to_2d=transform,
        plot_rgb_triangle=False,
        mesh_resolution=n,
    )
    # plt.xlim(-0.2, 0.9)
    # plt.ylim(+0.0, 0.7)
    plt.savefig(f"macadam-{n:03d}.png")
    return


def _parse_cmd_arguments():
    parser = argparse.ArgumentParser(
        description="Show piecewise linear transformation."
    )
    parser.add_argument("infile", help="input data file")
    return parser.parse_args()


if __name__ == "__main__":
    _main()
