#!/bin/bash
#TrialNum=20
#sleep 24h
model="'fmincon'"
#for N in $(seq 0.8 0.1 0.8) 
for TrialNum in $(seq 16 1 20)
do
for seperation in $(seq 0 1 7)
do
		     echo "#!/bin/bash">subjob.bash
		     echo '#SBATCH --job-name='$model,$TrialNum,$seperation>>subjob.bash
  	             echo '#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j' >> subjob.bash
   		     echo '#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j' >>subjob.bash		     
      		     echo "#SBATCH --exclusive">>subjob.bash
		     echo "#SBATCH --partition=fullnode">>subjob.bash
       		     echo "#SBATCH -N 1">>subjob.bash
		     echo "cd '/home/zhang.chi9/research/logscale/code/matlab'">>subjob.bash
       		     echo "matlab -nodisplay -r \"generate_network_N800($model,$TrialNum,$seperation)\"">>subjob.bash
       		     sbatch subjob.bash
	             sleep 30s
#sleep 12h
done
done
