#!/bin/bash
#TrialNum=20
#sleep 24h
N=800
for alpha in $(seq 0.8 0.1 0.8) 
do
		     echo "#!/bin/bash">subjob.bash
		     echo '#SBATCH --job-name=perceptron'$alpha,$N>>subjob.bash
  	             echo '#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j' >> subjob.bash
   		     echo '#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j' >>subjob.bash		     
      		     echo "#SBATCH --exclusive">>subjob.bash
		     echo '#SBATCH --partition=fullnode'>>subjob.bash
       		     echo "#SBATCH -N 1">>subjob.bash
		     #echo "#SBATCH --time=$( t2sd )">>subjob.bash
		     echo "cd '/home/zhang.chi9/research/logscale/code/matlab'">>subjob.bash
       		     echo 'matlab -nodisplay -r "perceptron_sub_prob('$N','$alpha')"'>>subjob.bash
       		     sbatch subjob.bash
	             #sleep 10s
#sleep 12h
done
