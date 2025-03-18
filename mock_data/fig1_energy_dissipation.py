import numpy as np
import matplotlib.pyplot as plt

# Parameters from Eq. (1)
t = np.linspace(0, 1e-12, 1000)  # Time (s)
E0 = 1e4  # Initial energy (eV)
Gamma_phys = 1e15  # High-energy damping
Gamma_chem = 1e12  # Chemical damping

E_phys = E0 * np.exp(-Gamma_phys * t)
E_chem = E0 * np.exp(-Gamma_chem * t)

plt.plot(t, E_phys, label='High-energy')
plt.plot(t, E_chem, label='Chemical')
plt.xlabel('Time (s)'); plt.ylabel('Energy (eV)')
plt.legend(); plt.savefig('Fig1_EnergyDissipation.png')
