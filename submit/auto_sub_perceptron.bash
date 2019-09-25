#!/bin/bash
#TrialNum=20
#sleep 24h
N=800
for Trialnum in $(seq 1 1 100) 
do
		     echo "#!/bin/bash">subjob.bash
		     echo '#SBATCH --job-name=perceptron'$Trialnum,$N>>subjob.bash
  	             echo '#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j' >> subjob.bash
   		     echo '#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j' >>subjob.bash		     
      		     echo "#SBATCH --exclusive">>subjob.bash
		     echo '#SBATCH --partition=fullnode'>>subjob.bash
       		     echo "#SBATCH -N 1">>subjob.bash
       		     echo 'matlab -nodisplay -r "perceptron_rule_generate_network('$N','$Trialnum')"'>>subjob.bash
       		     sbatch subjob.bash
	             #sleep 10s
#sleep 12h
done
