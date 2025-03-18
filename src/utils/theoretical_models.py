def energy_dissipation(t, E0=1e4, Gamma=1e12):
    """Eq.1: Exponential energy dissipation (mock dynamics)"""
    return E0 * np.exp(-Gamma * t)

def tfm_orbital_energy(n, l, lambda_beta_sq=1e-5):
    """Eq.2: TFM-corrected orbital energies"""
    E_QM = -13.6 / n**2
    f_nl = (1 + 0.1/n**2) + l*(l+1)*1e-3
    return E_QM * (1 + lambda_beta_sq * f_nl)

def reaction_rate(t, k_std=1e-3, Gamma_chem=1e12, A_osc=0.01):
    """Eq.6: Oscillatory reaction rates (mock coherence)"""
    return k_std * np.exp(-Gamma_chem*t) * (1 + A_osc*np.cos(1e12*t))
