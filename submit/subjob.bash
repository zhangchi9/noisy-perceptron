#!/bin/bash
#SBATCH --job-name=perceptron0.8,800
#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j
#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j
#SBATCH --exclusive
#SBATCH --partition=fullnode
#SBATCH -N 1
cd '/home/zhang.chi9/research/logscale/code/matlab'
matlab -nodisplay -r "perceptron_sub_prob(800,0.8)"
