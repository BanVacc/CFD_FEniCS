import matplotlib.pyplot as plt
import mshr as ms
from fenics import *
from solvers.SolverBase import SolverBase

from Fluid import Fluid
from Problem import Problem


class NewtonSteadySolver(SolverBase):
    def Solve(self,
              fluid: Fluid,
              problem: Problem,
              final_time: float,
              filename: str,
              time_step: float = 10 ** -3):

        V = VectorElement("P", problem.mesh.ufl_cell(), 2)
        Q = FiniteElement("P", problem.mesh.ufl_cell(), 1)
        TH = MixedElement([V, Q])
        W = FunctionSpace(problem.mesh, TH)

        w = Function(W)
        u, p = split(w)
        (v, q) = TestFunctions(W)

        f = Constant((0.0, 0.0))
        nu = Constant(fluid.mu/fluid.rho)
        F = (inner(grad(u)*u, v)*dx
             + nu*inner(grad(u), grad(v))*dx
             - inner(p, div(v))*dx
             + inner(q, div(u))*dx
             + inner(f, v)*dx)

        with XDMFFile(MPI.comm_world, f'{filename}.xdmf') as file:
            file.parameters.update(
                {
                    "functions_share_mesh": True,
                    "rewrite_function_mesh": False
                })

            for ic in problem.ics:  # todo
                w.interpolate(ic)

            bcs = []
            for bc in problem.bcs:
                subspace = bc[0]
                expression = bc[1]
                condition = bc[2]
                bcs.append(DirichletBC(W.sub(subspace), expression, condition))

            solve(F == 0, w, bcs,
                  solver_parameters={"newton_solver": {
                      "relative_tolerance": 1e-6,
                      "maximum_iterations": 15}})

            (u, p) = w.split(deepcopy=True)
            u.rename('U', 'Velocity')
            p.rename('P', 'Pressure')
            file.write(u, 0)
            file.write(p, 0)

# # Посмотреть вариант сдесь
# # https://fenicsproject.org/qa/13947/newton-solver-did-not-converge/
