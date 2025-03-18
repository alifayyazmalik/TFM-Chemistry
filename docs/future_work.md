# How to Run Real Simulations with TFM

## Step 1: Adapt Templates for Your HPC
- Use `templates/fenics_pde.py` to solve Eqs.12-14.  
- Replace mock functions in `src/theory/` with FEniCS/PETSc solvers.

## Step 2: Validate Against Experiments
1. **Atomic Spectroscopy**:  
   - Compare `src/theory/orbital_energy.py` predictions to Rydberg state data.  
   - Required precision: \( \Delta E/E \sim 10^{-5} \).  

2. **Ultra-Cold Reactions**:  
   - Test oscillatory rates (`src/theory/reaction_rates.py`) in molecular beam experiments.  

## Step 3: Contribute Back
- Share your HPC/experimental results via GitHub Issues or Pull Requests!
