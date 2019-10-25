#!/bin/bash
#TrialNum=20
#sleep 24h
model="'fmincon'"
#for N in $(seq 0.8 0.1 0.8) 
for N in 200 
do
for i in $(seq 1 1 10)
do
		     echo "#!/bin/bash">subjob.bash
		     echo '#SBATCH --job-name='$model,$N>>subjob.bash
  	             echo '#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j' >> subjob.bash
   		     echo '#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j' >>subjob.bash		     
      		     echo "#SBATCH --exclusive">>subjob.bash
		     echo "#SBATCH --partition=fullnode">>subjob.bash
       		     echo "#SBATCH -N 1">>subjob.bash
		     echo "cd '/home/zhang.chi9/research/logscale/code/functions'">>subjob.bash
       		     echo "matlab -nodisplay -r \"generate_network($N,-6,40,$model,'/scratch/zhang.chi9/perceptron/data/different_network_size_load_at_capacity')\"">>subjob.bash
       		     sbatch subjob.bash
	             sleep 30s
#sleep 12h
done
done
