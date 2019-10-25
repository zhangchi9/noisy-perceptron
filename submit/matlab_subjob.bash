#!/bin/bash
#SBATCH --job-name='fmincon',30,7
#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j
#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j
#SBATCH --exclusive
#SBATCH --partition=fullnode
#SBATCH -N 1
cd '/home/zhang.chi9/research/logscale/code/matlab'
matlab -nodisplay -r "generate_network_N800('fmincon',30,7)"
