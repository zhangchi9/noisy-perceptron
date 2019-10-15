#!/bin/bash
#TrialNum=20
#sleep 24h
N=400
for rj in $(seq -4 -0.5 -12) 
do
        for a in 5   10    15    20    25 
        do
		     echo "#!/bin/bash">subjob.bash
		     echo '#SBATCH --job-name=generate_nets'$rj,$a>>subjob.bash
  	             echo '#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j' >> subjob.bash
   		     echo '#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j' >>subjob.bash		     
      		     echo "#SBATCH --exclusive">>subjob.bash
		     echo '#SBATCH --partition=fullnode'>>subjob.bash
       		     echo "#SBATCH -N 1">>subjob.bash
		     #echo "#SBATCH --time=$( t2sd )">>subjob.bash
		     echo "cd '/home/zhang.chi9/research/logscale/code/matlab'">>subjob.bash
       		     echo 'matlab -nodisplay -r "fmincon_generate_network2('$rj','$a','$N')"'>>subjob.bash
       		     sbatch subjob.bash
	             sleep 10s
	done
#sleep 12h
done
