#!/bin/bash
#SBATCH --job-name=retrieval_no_noise
#SBATCH --output=out.%j
#SBATCH --error=err.%j
#SBATCH -N 1
#SBATCH -p fullnode
#SBATCH --exclusive
matlab -nodisplay -r prob_retrieval_map
