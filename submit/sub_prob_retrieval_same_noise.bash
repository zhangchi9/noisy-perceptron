#!/bin/bash
#SBATCH --job-name=retrieval_map
#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j
#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j
#SBATCH -N 1
#SBATCH -p fullnode
#SBATCH --exclusive
cd '/home/zhang.chi9/research/logscale/code/matlab'
matlab -nodisplay -r "prob_retrieval_map()"
