#!/bin/bash
#SBATCH --job-name=dynamics
#SBATCH --output=out.%j
#SBATCH --error=err.%j
#SBATCH --exclusive
#SBATCH -p fullnode
#SBATCH -N 1
matlab -nodisplay -r dynamics_analysis_add_noise
