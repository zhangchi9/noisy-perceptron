
function [zscore_mean,zscore_std] = get_zscore2(rin,a,trial_num_range)

zscore = nan(length(trial_num_range),2);
for i = 1:length(trial_num_range)
    trial_num = trial_num_range(i);
    fname = ['rin_',num2str(rin),'a_',num2str(a),'TrialNum_',num2str(trial_num),'.mat'];
    if ~isfile(fname)
        zscore(i,:) = nan(1,2);
    else
        load(fname)
        shuffled_time = 500;
        W = W';
        motif22E = nan(shuffled_time,2);
        %make_motif34lib
        N = 400;
        Ninh = 80;
        
        neuron_id = sum(abs(W(Ninh+1:end,:)),2)<N*71;%check if a correct solution
        WE = W(Ninh+1:end,Ninh+1:end);
        WE(WE<3) = 0; %
        WE(WE~=0)=1;
        WE = WE(neuron_id,neuron_id);
        motif21E = motif2struct_bin(WE);
        
        for n = 1:shuffled_time
            W_shuffle = dir_generate_srand(WE);
            motif22E(n,:) = motif2struct_bin(W_shuffle);
        end
        
        z = (motif21E' - mean(motif22E))./std(motif22E);
        znorm = z/sqrt(sum(z.^2));
        zscore(i,:) = znorm;
    end
end
zscore_mean = nanmean(zscore);
zscore_std = nanstd(zscore);
end

function count = motif2struct_bin(WE)
WW = WE.*WE';
count = zeros(2,1);
n_uc = sum(sum(WE))-sum(WW(:)); %number unidirectional connections
n_bc = (sum(WW(:)) - sum(diag(WW)))/2;%number of bidirectional connections
count(1,1) = n_uc;
count(2,1) = n_bc;
end

