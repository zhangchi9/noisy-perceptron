#!/bin/bash
#SBATCH --job-name=perceptron100,800
#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j
#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j
#SBATCH --exclusive
#SBATCH --partition=fullnode
#SBATCH -N 1
matlab -nodisplay -r "perceptron_rule_generate_network(800,100)"
