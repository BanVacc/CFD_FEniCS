from fenics import *
class Problem():
    def __init__(self,mesh,ics,bcs) -> None:
        self.mesh = mesh
        self.ics = ics
        self.bcs = bcs

        
            
