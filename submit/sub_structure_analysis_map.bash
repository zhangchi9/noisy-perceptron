#!/bin/bash
#SBATCH --job-name=structure_analysis
#SBATCH --output=out.%j
#SBATCH --error=err.%j
#SBATCH -N 1
#SBATCH -p fullnode
#SBATCH --exclusive
matlab -nodisplay -r "structure_analysis_map_new_fit()"
