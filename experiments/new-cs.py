# -*- coding: utf-8 -*-
#
from __future__ import print_function, division

import tempfile
import os

from dolfin import (
    Mesh, FunctionSpace, Function, grad, VectorFunctionSpace, project,
    TrialFunction, TestFunction, dot, dx, assemble
    )
import matplotlib.pyplot as plt
import numpy
from scipy.sparse.linalg import LinearOperator
from scipy.optimize import leastsq, least_squares, minimize
import yaml

import colorio
import meshio
import meshzoo

from pade2d import Pade2d


def f_ellipse(a_b_theta, x):
    a, b, theta = a_b_theta
    cos = numpy.cos(theta)
    sin = numpy.sin(theta)
    return (
        + a**2 * (x[0]*cos + x[1]*sin)**2
        + b**2 * (x[0]*sin - x[1]*cos)**2
        - 1.0
        )


def jac_ellipse(a_b_theta, x):
    a, b, theta = a_b_theta
    cos = numpy.cos(theta)
    sin = numpy.sin(theta)
    return numpy.array([
        + 2*a * (x[0]*cos + x[1]*sin)**2,
        #
        + 2*b * (x[0]*sin - x[1]*cos)**2,
        #
        + a**2 * 2*(x[0]*cos + x[1]*sin) * (-x[0]*sin + x[1]*cos)
        + b**2 * 2*(x[0]*sin - x[1]*cos) * (+x[0]*cos + x[1]*sin),
        ]).T


def _get_luo_rigg():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, '../colorio/data/luo-rigg/luo-rigg.yaml')) as f:
        data = yaml.safe_load(f)

    centers = []
    J = []
    for _, data_set in data.items():
        for _, dat in data_set.items():
            x, y, Y, a, ab, theta, _ = dat
            a /= 1.0e4
            a *= (Y/30)**0.2
            b = a / ab

            centers.append([x, y])

            J.append(numpy.array([
                [a * numpy.cos(theta), -b * numpy.sin(theta)],
                [a * numpy.sin(theta), b * numpy.cos(theta)],
                ]))

    return numpy.array(centers), numpy.moveaxis(numpy.array(J), 0, -1)


def _get_macadam():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, '../colorio/data/macadam1942/table3.yaml')) as f:
        data = yaml.safe_load(f)

    centers = []
    points = []
    for datak in data:
        # collect ellipse points
        _, _, _, _, delta_y_delta_x, delta_s = numpy.array(datak['data']).T
        if len(delta_s) < 2:
            continue
        center = [datak['x'], datak['y']]
        centers.append(center)
        offset = (
            numpy.array([numpy.ones(delta_y_delta_x.shape[0]), delta_y_delta_x])
            / numpy.sqrt(1 + delta_y_delta_x**2) * delta_s
            )
        points.append(numpy.column_stack([
            (center + offset.T).T,
            (center - offset.T).T,
            ]))

    centers = numpy.array(centers)
    J = get_local_linearizations1(centers, points)
    return centers, numpy.moveaxis(J, 0, -1)
    # return centers, self.get_local_linearizations2(centers, points)


def get_local_linearizations1(centers, points):
    # Get ellipse parameters
    X = [
        (pts.T - center).T
        for center, pts in zip(centers, points)
        ]
    a_b_theta = numpy.array([
        # Solve least squares problem for [1/a, 1/b, theta]
        # and pick [a, b, theta]
        leastsq(
            lambda a_b_theta: f_ellipse(a_b_theta, x),
            [1.0, 1.0, 0.0],
            Dfun=lambda a_b_theta: jac_ellipse(a_b_theta, x),
            )[0]
        for x in X
        ])
    a_b_theta = numpy.array([
        1 / a_b_theta[:, 0],
        1 / a_b_theta[:, 1],
        a_b_theta[:, 2]
        ]).T
    # Construct 2x2 matrices that approximately convert unit circles into
    # the ellipse defined by the points.
    J = []
    for abt in a_b_theta:
        a, b, theta = abt
        J.append(numpy.array([
            [a * numpy.cos(theta), -b * numpy.sin(theta)],
            [a * numpy.sin(theta), b * numpy.cos(theta)],
            ]))

    return numpy.array(J)


def get_local_linearizations2(centers, points):
    X = [
        (pts.T - center).T
        for center, pts in zip(centers, points)
        ]

    def f_linear_function(j, x):
        Jx = numpy.dot(j.reshape(2, 2), x)
        out = numpy.einsum('ij,ij->j', Jx, Jx) - 1.0
        return out

    def jac_linear_function(j, x):
        J = j.reshape(2, 2)
        return numpy.array([
            2*J[0, 0]*x[0]**2 + 2*J[0, 1]*x[0]*x[1],
            2*J[0, 1]*x[1]**2 + 2*J[0, 0]*x[0]*x[1],
            2*J[1, 0]*x[0]**2 + 2*J[1, 1]*x[0]*x[1],
            2*J[1, 1]*x[1]**2 + 2*J[1, 0]*x[0]*x[1],
            ]).T

    J = []
    for x in X:
        j, _ = leastsq(
            lambda J: f_linear_function(J, x),
            [1.0, 0.0, 0.0, 1.0],
            Dfun=lambda J: jac_linear_function(J, x),
            # full_output=True
            )
        J.append(numpy.linalg.inv(j.reshape(2, 2)))

    return numpy.array(J)


class PadeEllipse(object):
    def __init__(self, centers, J, degrees):
        self.centers = centers
        self.J = J

        self.target = 0.002
        self.J /= self.target

        self.num_f_eval = 0

        self.degrees = degrees

        num_coefficients = [
            (degrees[0]+1) * (degrees[0]+2) // 2,
            (degrees[1]+1) * (degrees[1]+2) // 2,
            (degrees[2]+1) * (degrees[2]+2) // 2,
            (degrees[3]+1) * (degrees[3]+2) // 2,
            ]

        # Choose the coefficiens to create the identity function
        ax = numpy.zeros(num_coefficients[0])
        ax[1] = 1
        bx = numpy.zeros(num_coefficients[1] - 1)
        ay = numpy.zeros(num_coefficients[2])
        ay[2] = 1
        by = numpy.zeros(num_coefficients[3] - 1)

        self.alpha = numpy.concatenate([ax, bx, ay, by])

        bx = numpy.concatenate([[1.0], bx])
        by = numpy.concatenate([[1.0], by])

        self.pade2d = Pade2d(self.centers.T, degrees, ax, bx, ay, by)

        # self.J = numpy.array(self.get_local_linearizations2(centers, points))

        # # plot
        # for center, pts, j in zip(centers, points, self.J):
        #     # plot points
        #     p = (pts.T - center).T
        #     plt.plot(*p, '.')
        #     # plot circle
        #     t = numpy.linspace(0.0, 2.0*numpy.pi, 1000)
        #     xy = numpy.array([numpy.cos(t), numpy.sin(t)])
        #     plt.plot(*numpy.dot(j, xy), '-', label='ellipse')
        #     plt.legend()
        #     # # plot transformation
        #     # xy_new = numpy.dot(j, p)
        #     # plt.plot(*xy_new, 'x')
        #     plt.axis('equal')
        #     plt.show()
        return

    def _set_alpha(self, alpha):
        # Subtract 1 for each denominator polynomial since the constant
        # coefficient is fixed to 1.0.
        assert len(alpha) == len(self.alpha)

        self.alpha = alpha

        num_coefficients = [(d+1)*(d+2)//2 for d in self.degrees]
        num_coefficients[1] -= 1
        num_coefficients[3] -= 1

        ax, bx, ay, by = \
            numpy.split(alpha, numpy.cumsum(num_coefficients[:-1]))
        bx = numpy.concatenate([[1.0], bx])
        by = numpy.concatenate([[1.0], by])

        self.pade2d.set_coefficients(ax, bx, ay, by)
        return

    def get_q2_r2(self, alpha):
        self._set_alpha(alpha)

        # jacs and J are of shape (2, 2, k). M must be of the same shape and
        # contain the result of the k 2x2 dot products. Perhaps there's a
        # dot() for this.
        M = numpy.einsum('ijl,jkl->ikl', self.pade2d.jac(), self.J)

        # One could use
        #
        #     M = numpy.moveaxis(M, -1, 0)
        #     _, sigma, _ = numpy.linalg.svd(M)
        #
        # but computing the singular values explicitly via
        # <https://scicomp.stackexchange.com/a/14103/3980> is faster.
        a = (M[0, 0] + M[1, 1]) / 2
        b = (M[0, 0] - M[1, 1]) / 2
        c = (M[1, 0] + M[0, 1]) / 2
        d = (M[1, 0] - M[0, 1]) / 2

        # From the square roots of q2 and r2, the ellipse axes can be computed,
        # namely
        #
        #   s1 = q + r
        #   s2 = q - r
        #
        q2 = a**2 + d**2
        r2 = b**2 + c**2

        return q2, r2

    def get_ellipse_axes(self, alpha):
        q, r = numpy.sqrt(self.get_q2_r2(alpha))
        sigma = numpy.array([q+r, q-r]) * self.target
        return sigma

    def cost(self, alpha):
        q2, r2 = self.get_q2_r2(alpha)

        out = numpy.array([q2 - 1.0, r2]).flatten()

        if self.num_f_eval % 10000 == 0:
            cost = numpy.sum(out**2)
            print('{:7d}     {}'.format(self.num_f_eval, cost))

        self.num_f_eval += 1
        return out


class PiecewiseEllipse(object):
    def __init__(self, centers, J):
        self.centers = centers
        self.J = J

        self.target = 0.002
        self.J /= self.target

        # dir_path = os.path.dirname(os.path.realpath(__file__))
        # with open(os.path.join(dir_path, '../colorio/data/gamut_triangulation.yaml')) as f:
        #     data = yaml.safe_load(f)

        # self.points = numpy.column_stack([
        #     data['points'], numpy.zeros(len(data['points']))
        #     ])
        # self.cells = numpy.array(data['cells'])

        # self.points, self.cells = colorio.xy_gamut_mesh(0.15)

        self.points, self.cells = meshzoo.triangle(
            corners=numpy.array([
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]
                ]),
            ref_steps=3
            )

        # https://bitbucket.org/fenics-project/dolfin/issues/845/initialize-mesh-from-vertices
        with tempfile.TemporaryDirectory() as temp_dir:
            tmp_filename = os.path.join(temp_dir, 'test.xml')
            meshio.write(
                tmp_filename, self.points, {'triangle': self.cells},
                file_format='dolfin-xml'
                )
            mesh = Mesh(tmp_filename)

        self.V = FunctionSpace(mesh, 'CG', 1)
        self.Vgrad = VectorFunctionSpace(mesh, 'DG', 0)

        self.ux = Function(self.V)
        self.uy = Function(self.V)

        # self.alpha = self.u.vector().get_local()
        self.ax = self.ux.vector().get_local()
        self.ay = self.uy.vector().get_local()
        self.alpha = numpy.concatenate([self.ax, self.ay])

        self.num_f_eval = 0

        u = TrialFunction(self.V)
        v = TestFunction(self.V)
        self.L = assemble(dot(grad(u), grad(v)) * dx)

        # The functions are locally linear, so we can cast the projection into a
        # matrix. The matrix will be very sparse, but never mind that now; just
        # use an array.
        n = self.V.dim()
        self.jacs = []
        for k in range(n):
            e = numpy.zeros(n)
            e[k] = 1.0

            u = Function(self.V)
            u.vector().set_local(e)

            j = project(grad(u), self.Vgrad)
            self.jacs.append([j(x, y) for x, y in self.centers])

        self.jacs = numpy.array(self.jacs).T
        return

    def get_q2_r2(self):
        # Keep an eye on
        # https://www.allanswered.com/post/eqkmm/matrix-from-project-evaluation-at-many-points-at-once/
        # https://bitbucket.org/fenics-project/dolfin/issues/1011/evaluate-expressions-at-many-points-at
        # for speeding this up.
        # jx = project(grad(self.ux), self.Vgrad)
        # jy = project(grad(self.uy), self.Vgrad)
        # jac = numpy.array([[jx(x, y), jy(x, y)] for x, y in self.centers])
        # jac = numpy.moveaxis(jac, 0, -1)
        # assert numpy.all(abs(jj - jac) < 1.0e-15)

        # u = numpy.array([
        #     self.ux.vector().get_local(),
        #     self.uy.vector().get_local(),
        #     ])
        # jac = numpy.einsum('ikl,jl->jik', self.jacs, u)

        jac_x = numpy.dot(self.jacs, self.ux.vector().get_local())
        jac_y = numpy.dot(self.jacs, self.uy.vector().get_local())
        jac = numpy.array([jac_x, jac_y])

        # jacs and J are of shape (2, 2, k). M must be of the same shape and
        # contain the result of the k 2x2 dot products. Perhaps there's a
        # dot() for this.
        M = numpy.einsum('ijl,jkl->ikl', jac, self.J)
        # M = numpy.array([
        #     [
        #         jac[0][0]*self.J[0][0] + jac[0][1]*self.J[1][0],
        #         jac[0][0]*self.J[0][1] + jac[0][1]*self.J[1][1],
        #     ],
        #     [
        #         jac[1][0]*self.J[0][0] + jac[1][1]*self.J[1][0],
        #         jac[1][0]*self.J[0][1] + jac[1][1]*self.J[1][1],
        #     ],
        #     ])

        # One could use
        #
        #     M = numpy.moveaxis(M, -1, 0)
        #     _, sigma, _ = numpy.linalg.svd(M)
        #
        # but computing the singular values explicitly via
        # <https://scicomp.stackexchange.com/a/14103/3980> is faster.
        a = (M[0, 0] + M[1, 1]) / 2
        b = (M[0, 0] - M[1, 1]) / 2
        c = (M[1, 0] + M[0, 1]) / 2
        d = (M[1, 0] - M[0, 1]) / 2

        # From the square roots of q2 and r2, the ellipse axes can be computed,
        # namely
        #
        #   s1 = q + r
        #   s2 = q - r
        #
        q2 = a**2 + d**2
        r2 = b**2 + c**2

        return q2, r2

    def get_ellipse_axes(self, alpha):
        ax, ay = numpy.split(alpha, 2)
        self.ux.vector().set_local(ax)
        self.uy.vector().set_local(ay)

        q, r = numpy.sqrt(self.get_q2_r2())
        sigma = numpy.array([q+r, q-r]) * self.target
        return sigma

    def cost(self, alpha):
        ax, ay = numpy.split(alpha, 2)
        self.ux.vector().set_local(ax)
        self.uy.vector().set_local(ay)

        res_x = self.L * self.ux.vector()
        res_y = self.L * self.uy.vector()

        q2, r2 = self.get_q2_r2()

        # scale the problem with 100 to avoid premature exit in optimization
        out = 100 * numpy.array([
            res_x.get_local() / len(res_x.get_local()),
            res_y.get_local() / len(res_y.get_local()),
            (q2 - 1.0) / len(q2),
            r2 / len(r2),
            ])

        if self.num_f_eval % 1000 == 0:
            cost = numpy.array([numpy.sum(ot**2) for ot in out])
            print('{:7d}     {:e} {:e} {:e} {:e}'.format(self.num_f_eval, *cost))

        self.num_f_eval += 1
        return numpy.concatenate(out)

    def get_jac(self, alpha):
        m = 2*self.V.dim() + 2*self.centers.shape[0]
        n = alpha.shape[0]

        ax, ay = numpy.split(alpha, 2)
        jac_alpha_x = numpy.dot(self.jacs, ax)
        jac_alpha_y = numpy.dot(self.jacs, ay)
        jac_alpha = numpy.array([jac_alpha_x, jac_alpha_y])
        M_alpha = numpy.einsum('ijl,jkl->ikl', jac_alpha, self.J)
        #
        a_alpha = (M_alpha[0, 0] + M_alpha[1, 1]) / 2
        b_alpha = (M_alpha[0, 0] - M_alpha[1, 1]) / 2
        c_alpha = (M_alpha[1, 0] + M_alpha[0, 1]) / 2
        d_alpha = (M_alpha[1, 0] - M_alpha[0, 1]) / 2

        # build lower q, r part of the linear operator
        dq2 = []
        dr2 = []
        for k in range(n):
            e = numpy.zeros(n)
            e[k] = 1.0

            ax, ay = numpy.split(e, 2)

            # q2, r2 part
            jac_phix = numpy.dot(self.jacs, ax)
            jac_phiy = numpy.dot(self.jacs, ay)
            jac_phi = numpy.array([jac_phix, jac_phiy])
            M_phi = numpy.einsum('ijl,jkl->ikl', jac_phi, self.J)
            #
            a_phi = (M_phi[0, 0] + M_phi[1, 1]) / 2
            b_phi = (M_phi[0, 0] - M_phi[1, 1]) / 2
            c_phi = (M_phi[1, 0] + M_phi[0, 1]) / 2
            d_phi = (M_phi[1, 0] - M_phi[0, 1]) / 2

            dq2.append((a_alpha*a_phi + d_alpha*d_phi) / 2)
            dr2.append((b_alpha*b_phi + c_alpha*c_phi) / 2)
        dq2 = numpy.column_stack(dq2)
        dr2 = numpy.column_stack(dr2)

        def matvec(phi):
            if len(phi.shape) > 1:
                assert len(phi.shape) == 2
                assert phi.shape[1] == 1
                phi = phi[:, 0]

            ax, ay = numpy.split(phi, 2)

            # Laplace part (it's linear, so this is easy)
            self.ux.vector().set_local(ax)
            self.uy.vector().set_local(ay)
            res_x = self.L * self.ux.vector()
            res_y = self.L * self.uy.vector()

            # q2, r2 part
            dq2_phi = numpy.dot(dq2, phi)
            dr2_phi = numpy.dot(dr2, phi)

            # scale the problem with 100 to avoid premature exit in optimization
            out = 100 * numpy.array([
                res_x.get_local() / len(res_x.get_local()),
                res_y.get_local() / len(res_y.get_local()),
                dq2_phi / len(dq2),
                dr2_phi / len(dr2),
                ])
            return numpy.concatenate(out)

        def rmatvec(vec):
            res_x, res_y, dq2_phi, dr2_phi = numpy.split(
                vec,
                numpy.cumsum([
                    self.V.dim(), self.V.dim(), self.centers.shape[0]
                    ])
                )
            res_x *= 100 / len(res_x)
            res_y *= 100 / len(res_x)
            dq2_phi *= 100 / len(dq2_phi)
            dr2_phi *= 100 / len(dr2_phi)

            self.ux.vector().set_local(res_x)
            self.uy.vector().set_local(res_y)
            phi = numpy.concatenate([
                (self.L * self.ux.vector()).get_local(),
                (self.L * self.uy.vector()).get_local(),
                ])

            q2p = numpy.dot(dq2.T, dq2_phi)
            r2p = numpy.dot(dr2.T, dr2_phi)
            return phi + q2p + r2p

        return LinearOperator([m, n], matvec=matvec, rmatvec=rmatvec)


def _main():
    centers, J = _get_macadam()
    # centers, J = _get_luo_rigg()

    # problem = PadeEllipse(centers, J, [2, 0, 2, 0])
    problem = PiecewiseEllipse(centers, J)

    print('num parameters: {}'.format(len(problem.alpha)))

    alpha0 = problem.alpha.copy()

    numpy.random.seed(0)
    alpha0 = numpy.random.rand(*alpha0.shape)

    # Levenberg-Marquardt (lm) is better suited for small, dense, unconstrained
    # problems, but it needs more conditions than parameters. This is not the
    # case for larger polynomial degrees.
    print('f evals     cost')
    out = least_squares(
        problem.cost, alpha0,
        jac=problem.get_jac,
        method='trf'
        )
    # out = minimize(problem.cost, alpha0, method='Powell')

    final_cost = numpy.sum(problem.cost(out.x)**2)
    print('{:7d}     {}'.format(problem.num_f_eval, final_cost))

    # # plot statistics
    # axes0 = problem.get_ellipse_axes(alpha0).T.flatten()
    # plt.plot(axes0, label='axes lengths before')
    # axes1 = problem.get_ellipse_axes(out.x).T.flatten()
    # plt.plot(axes1, label='axes lengths opt')
    # plt.legend()
    # plt.grid()

    # # Plot unperturbed MacAdam
    # plt.figure()
    # # colorio.plot_luo_rigg(
    # #     ellipse_scaling=1,
    # colorio.plot_macadam(
    #     ellipse_scaling=10,
    #     plot_rgb_triangle=False,
    #     )

    # Plot perturbed MacAdam
    def transform(XY):
        is_solo = len(XY.shape) == 1
        if is_solo:
            XY = numpy.array([XY]).T
        # print(XY)
        out = numpy.array([
            [problem.ux(x, y) for x, y in XY.T],
            [problem.uy(x, y) for x, y in XY.T],
            ])
        if is_solo:
            out = out[..., 0]
        return out

    plt.figure()
    # colorio.plot_luo_rigg(
    #     ellipse_scaling=1,
    colorio.plot_macadam(
        ellipse_scaling=10,
        # xy_to_2d=problem.pade2d.eval,
        xy_to_2d=transform,
        plot_rgb_triangle=False,
        )

    plt.show()
    return


if __name__ == '__main__':
    _main()
