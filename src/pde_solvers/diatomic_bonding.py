from fenics import *
import matplotlib.pyplot as plt

def diatomic_bond_energy(r, lambda_beta_sq=1e-5):
    """
    Morse-like bonding energy (Eq. 4).
    """
    return - (1/r) * (1 - np.exp(-lambda_beta_sq * r))

# Generate mock data for Fig2
r = np.linspace(0.5, 5, 100)  # Bond length (Å)
E_bond = [diatomic_bond_energy(ri) for ri in r]
