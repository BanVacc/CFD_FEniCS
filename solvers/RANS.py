# Reynolds Average Navier-Stoke

# Kolmagorov scales
# length scale: eta = theta*(Re^(-3/4))
# time scale: tau_eta = theta* Re^(-1/2)
# velocity: u_eta = theta * (Re^(-1/4))
# 
# In order to converge
# dx < eta | dx 
# dt < tau_eta
# ?

from solvers.SolverBase import SolverBase
from Fluid import Fluid
from Problem import Problem
from fenics import *

class RANSSolver(SolverBase):

    def Solve(self,
               fluid: Fluid, 
               problem: Problem, 
               final_time: float, 
               filename: str, 
               time_step: float = 10 ** -3):
        
        Re=1
        
        # Compute Kolmagorov scales
        eta = Re**(-3/4)
        tau_eta = Re**(-1/2)
        u_eta = Re**(-1/4)

        print('Required number of steps:',1/tau_eta) # ? 21:00 YouTube:Why are Direct Numerical Simulations often impossible?

        if CellDiameter(problem.mesh) > eta:
            print('Will not converge')
            return
        if time_step > tau_eta:
            print('Will not converge')
            return