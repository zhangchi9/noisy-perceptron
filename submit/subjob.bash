#!/bin/bash
#SBATCH --job-name=perceptron0.40,200
#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j
#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j
#SBATCH --exclusive
#SBATCH --partition=fullnode
#SBATCH -N 1
cd '/home/zhang.chi9/research/logscale/code/matlab'
matlab -nodisplay -r "perceptron_sub_prob(200,0.40)"
