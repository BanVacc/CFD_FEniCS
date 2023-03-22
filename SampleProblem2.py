from Problem import Problem
import mshr as ms
from fenics import *
from Fluid import Fluid
from solvers.Chorin import ChorinSolver

parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"

mesh = UnitSquareMesh(41, 41, 'crossed')

ics = []
walls = 'on_boundary'
top_boundary = 'on_boundary && near(x[1],1.0)'
bottom_boundary = 'on_boundary && near(x[1],0.0)'
left_boundary = 'on_boundary && near(x[0],0.0)'
right_boundary = 'on_boundary && near(x[0],1.0)'
bcs = [
    [0, Constant((0.0, 0.0)), walls],  # No-slip
    [0, Constant((1.0, 0.0)), top_boundary],  # Velocity's x = 1 at top
    [1, Constant(0.0), top_boundary],  # Pressure out at top
    [2, Constant(15), left_boundary],
    [2, Constant(0), right_boundary]
]

problem = Problem(mesh, ics, bcs)
fluid = Fluid(1.0, 0.01, 0.598, 4184)

steady_solver = ChorinSolver()
steady_solver.Solve(fluid, problem, 1, 'LidDrivenCavityChorin', 0.01)
