from Problem import Problem
import mshr as ms
from fenics import *
from Fluid import Fluid
from solvers.NewtonSteady import NewtonSteadySolver


mesh = ms.generate_mesh(
    ms.Rectangle(
        Point(0.0, 0.0),
        Point(2.2, 0.41))
    -
    ms.Circle(Point(0.2, 0.2), 0.05), 96*2)

ics = []
bcs = [
    [0, Constant((0.0, 0.0)), 'on_boundary'],
    [0, Expression(
        ('4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)', '0'),
        degree=2), 'on_boundary && near(x[0], 0.0)'],
    [1, Constant(0.0), 'on_boundary && near(x[0], 2.2)']]

problem = Problem(mesh, ics, bcs)
fluid = Fluid(1.0, 0.01, 1, 1)

steady_solver = NewtonSteadySolver()
steady_solver.Solve(fluid, problem, 1, 'fac', 10**-3)
