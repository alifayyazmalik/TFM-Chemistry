import numpy as np
from src.utils import theoretical_models as tm

# Mock 3-atom system with random distances
positions = np.random.rand(3, 3)  # Fake coordinates
coherence = 0.05 * np.exp(-1.0 * np.linalg.norm(positions, axis=1))

# Plotting "stability map" (no actual PDE solved)
plt.scatter(positions[:,0], positions[:,1], c=coherence)
plt.savefig('Fig5_MultiAtom.png')
