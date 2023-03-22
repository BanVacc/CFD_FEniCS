from abc import ABC, abstractmethod

from Fluid import Fluid
from Problem import Problem


class SolverBase(ABC):

    @abstractmethod
    def Solve(self,
              fluid: Fluid,
              problem: Problem,
              final_time: float,
              filename:str,
              time_step: float = 10**-3
              ):
        pass
