#!/bin/bash
while true; do
 #sbatch sub_dynamics_analysis_add_noise.bash
 sbatch sub_structure_analysis_map.bash
 sbatch sub_prob_retrieval_same_noise.bash
 sleep 24h

done
