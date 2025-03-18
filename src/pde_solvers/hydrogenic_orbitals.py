import numpy as np
from fenics import *

def compute_tfm_energy(n, l, lambda_beta_sq=1e-5):
    """
    Compute TFM-corrected hydrogenic energy levels (Eq. 2).
    """
    E_QM = -13.6 / n**2  # Base QM energy (eV)
    f_nl = (1 + 0.1/n**2) + l*(l+1)*1e-3
    E_TFM = E_QM * (1 + lambda_beta_sq * f_nl)
    return E_TFM

# Example usage for high-n Rydberg states
n_values = np.arange(10, 100)
E_tfm = [compute_tfm_energy(n, l=0) for n in n_values]
