
function find_missing_submission(n_network)

file_obj = dir('/scratch/zhang.chi9/perceptron/data/tmp_N_800');

filename_list = {file_obj.name};

model = {'perceptron','fmincon'};

for mod = model
    for i = 1:n_network
        for j = 0:7
            filename = [mod{1},'_',num2str(i),'_',num2str(j),'.mat'];
            if ~any(strcmp(filename_list,filename))
                submit_sbatch(mod{1},i,j)
            end
        end
    end
end

function submit_sbatch(model,i,j)

fid = fopen('/home/zhang.chi9/research/logscale/code/submit/matlab_subjob.bash','w');

fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'#SBATCH --job-name=''%s'',%d,%d\n',model,i,j);
fprintf(fid,'#SBATCH --output=/home/zhang.chi9/research/logscale/logs/out.%%j\n');
fprintf(fid,'#SBATCH --error=/home/zhang.chi9/research/logscale/logs/err.%%j\n');
fprintf(fid,'#SBATCH --exclusive\n');
fprintf(fid,'#SBATCH --partition=fullnode\n');
fprintf(fid,'#SBATCH -N 1\n');
fprintf(fid,"cd '/home/zhang.chi9/research/logscale/code/matlab'\n");
fprintf(fid,'matlab -nodisplay -r "generate_network_N800(''%s'',%d,%d)"\n',model,i,j);

system('sbatch /home/zhang.chi9/research/logscale/code/submit/matlab_subjob.bash')













