import mshr as ms
from fenics import *

from Fluid import Fluid
from Problem import Problem
from solvers.Chorin import ChorinSolver
from solvers.NewtonSteady import NewtonSteadySolver

domain = ms.generate_mesh(
    ms.Rectangle(
        Point(0.0, 0.0),
        Point(2.2, 0.41))
    -
    ms.Circle(Point(0.2, 0.2), 0.05), 35)

ics = []
bcs = [
    [0, Constant((0.0, 0.0)), 'on_boundary'],
    [0, Expression(
        ('4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)', '0'),
        degree=2), 'on_boundary && near(x[0], 0.0)'],
    [1, Constant(0.0), 'on_boundary && near(x[0], 2.2)']]

problem = Problem(domain, ics, bcs)
fluid = Fluid(1.0, 0.001, 1, 1)

steady_solver = ChorinSolver()
steady_solver.Solve(fluid, problem, 1, 'Flow', 10**-3)
