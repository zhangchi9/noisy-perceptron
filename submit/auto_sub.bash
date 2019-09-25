#!/bin/bash
#TrialNum=20
#sleep 24h
for rj in $(seq -9.5 -0.5 -12) 
do
        for a in 0.01 5    10    15    20    25    30    35    40    45    50   55    60    65    70    75    80    85    90    95   100
        do
		     echo "#!/bin/bash">subjob.bash
		     echo '#SBATCH --job-name=generate_nets'$rj,$a>>subjob.bash
  	             echo '#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%j' >> subjob.bash
   		     echo '#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%j' >>subjob.bash		     
      		     echo "#SBATCH --exclusive">>subjob.bash
		     echo '#SBATCH --partition=fullnode'>>subjob.bash
       		     echo "#SBATCH -N 1">>subjob.bash
		     echo "#SBATCH --time=48:00:00" >>subjob.bash
       		     echo 'matlab -nodisplay -r "fmincon_generate_network2('$rj','$a')"'>>subjob.bash
       		     #sbatch subjob.bash
	             #sleep 10s
	done
#sleep 12h
done
