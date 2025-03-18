from fenics import *
import numpy as np

class MultiAtomSolver:
    def __init__(self, N_atoms, C=0.05, alpha=1.0):
        self.N = N_atoms
        self.C = C  # Coherence strength (eV)
        self.alpha = alpha  # Range parameter

    def wave_lump_pde(self, mesh):
        """
        Solve TFM PDE for N-atom system (Eq. 12-14).
        """
        V = FunctionSpace(mesh, 'P', 1)
        T_plus = TrialFunction(V)
        v = TestFunction(V)
        
        # Simplified quasi-static approximation (Eq. 14)
        a = dot(grad(T_plus), grad(v)) * dx
        L = self.C * exp(-self.alpha * T_plus) * v * dx
        T_sol = Function(V)
        solve(a == L, T_sol)
        return T_sol
