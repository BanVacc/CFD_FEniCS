from fenics import *
from Problem import Problem
from solvers.SolverBase import SolverBase
from tqdm import tqdm

from Fluid import Fluid


class ChorinSolver(SolverBase):
    def __init__(self) -> None:
        pass

    def Solve(self,
              fluid: Fluid,
              problem: Problem,
              final_time: float,
              filename: str,
              time_step: float = 10 ** -3):

        # Define function spaces
        V = VectorFunctionSpace(problem.mesh,'Lagrange',2)
        Q = FunctionSpace(problem.mesh,'Lagrange',1)

        # Define the trial and test functions
        u = TrialFunction(V)
        v = TestFunction(V)
        p = TrialFunction(Q)
        q = TestFunction(Q)
        t = TrialFunction(Q)
        s = TestFunction(Q)

        # todo: Define initial conditions
        pass

        # Define boundary conditions
        ubc = []
        pbc = []
        tbc = []

        for bc in problem.bcs:
            field = bc[0]
            expression = bc[1]
            condition = bc[2]

            if field == 0:
                ubc.append(DirichletBC(V, expression, condition))
            elif field == 1:
                pbc.append(DirichletBC(Q, expression, condition))
            elif field == 2:
                tbc.append(DirichletBC(Q, expression, condition))

        # Define the solution fields involved
        u_prev = Function(V) # velocity at previous time step, used for time integration
        u_tent = Function(V) # velocity before pressure correction
        u_next = Function(V,name='U') # velocity after pressure correction, at current time step
        p_next = Function(Q,name='P') # pressure at current time step
        t_prev = Function(Q)
        t_next = Function(Q,name='T')

        # Define constants
        dt = Constant(time_step)
        rho = Constant(fluid.rho)
        nu = Constant(fluid.mu / fluid.rho)
        kappa = Constant(fluid.alpha/(fluid.rho*fluid.cp))

        # Weak form of the momentum equation
        F1 = (
            (1.0 / dt) * inner(u - u_prev, v) * dx
            + inner(grad(u_prev) * u_prev, v) * dx  # Advection term
            + nu * inner(grad(u), grad(v)) * dx)  # Diffusion term

        # Weak form of the pressure poisson problem
        F2 = (
            inner(grad(p), grad(q)) * dx
            + (rho / dt) * div(u_tent) * q * dx)

        # Weak form of the velocity update equation
        F3 = (
            inner(u, v) * dx
            - inner(u_tent, v) * dx
            + (dt/rho) * inner(grad(p_next), v) * dx)

        # Weak form of heat transfer equation
        F4 = (
            (1.0 / dt) * inner(t - t_prev, s) * dx
            + dot(u_next, grad(t_prev)) * s * dx  # Advection
            + kappa * inner(grad(t_prev), grad(s)) * dx  # Diffusion term
        )

        # Split equations into lhs and rhs
        a1, L1 = system(F1)
        a2, L2 = system(F2)
        a3, L3 = system(F3)
        a4, L4 = system(F4)
        
        # Assemble matricies
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)
        A4 = assemble(a4)

        # Applying BCs to A matrecies
        [bc.apply(A1) for bc in ubc]
        [bc.apply(A2) for bc in pbc]
        [bc.apply(A3) for bc in ubc]
        [bc.apply(A4) for bc in tbc]

        # Time loop
        # Creating file for saving results
        with XDMFFile(MPI.comm_world, f'{filename}.xdmf') as file:
            file.parameters.update(
                {
                    "functions_share_mesh": True,
                    "rewrite_function_mesh": False
                })
            # Saving initial states
            file.write(u_next,0)
            file.write(p_next,0)
            file.write(t_next,0)

            # Use amg preconditioner if available
            prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

            n_iters = final_time // time_step
            # Simulation
            for i in tqdm(range(int(n_iters))):
                # (1) Solve for tentative velocity
                b1 = assemble(L1)
                [bc.apply(b1) for bc in ubc]
                solve(A1, u_tent.vector(), b1)

                # (2) Solve for the pressure
                b2 = assemble(L2)
                [bc.apply(b2) for bc in pbc]
                solve(A2, p_next.vector(), b2)

                # (3) Correct the velocities to be incompressible
                b3 = assemble(L3)
                [bc.apply(b3) for bc in ubc]
                solve(A3, u_next.vector(), b3)

                # (4) Solve heat advection-diffusion equation
                b4 = assemble(L4)
                [bc.apply(b4) for bc in tbc]
                solve(A4, t_next.vector(), b4)

                # Save to file
                file.write(u_next, i*time_step)
                file.write(p_next, i*time_step)
                file.write(t_next, i*time_step)

                # Advance in time
                t_prev.assign(t_next)
                u_prev.assign(u_next)
