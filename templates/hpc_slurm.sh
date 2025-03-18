#!/bin/bash
#SBATCH --job-name=TFM_HPC      # Rename for your simulation
#SBATCH --nodes=4               # ★ Modify for your cluster
#SBATCH --ntasks-per-node=128   # ★ Adjust based on core count

# ★★★ Replace with your actual PDE solver ★★★
module load fenics
mpirun -np 512 python your_actual_solver.py \
  --gamma_chem 1e12 \
  --lambda_beta_sq 1e-5 \
  --output real_simulation.h5
