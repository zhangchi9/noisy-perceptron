function structure_analysis_map(binwidth)
rin_ind_range = [-2:-0.5:-12];
a_range = [0.1,5:5:100];
trial_range = [1:100];
Pcon_inh = nan(length(rin_ind_range),length(a_range));
Pcon_exc = nan(length(rin_ind_range),length(a_range));
CV_inh = nan(length(rin_ind_range),length(a_range));
CV_exc = nan(length(rin_ind_range),length(a_range));
    
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
scratch_dir = ['/home/zhang.chi9/matlabtmp/', num2str(randi(10^9))];
mkdir(scratch_dir);
pc = parcluster('local');
pc.JobStorageLocation = scratch_dir;
parpool(pc,pc.NumWorkers);
parfor i = 1:length(rin_ind_range)
    tmp1 = nan(1,length(a_range));
    tmp2 = nan(1,length(a_range));
    tmp3 = nan(1,length(a_range));
    tmp4 = nan(1,length(a_range));
    for j = 1:length(a_range)
        [rin_ind_range(i),a_range(j)]
        [tmp1(j),tmp2(j),tmp3(j),tmp4(j)] = connection_prob_CV_averaged(rin_ind_range(i),a_range(j),trial_range,binwidth);
    end
    Pcon_inh(i,:) = tmp1;
    Pcon_exc(i,:) = tmp2;
    CV_inh(i,:) = tmp3;
    CV_exc(i,:) = tmp4;
end
save(['fmincon_structure_map_load_at_numerical_capacity','binwidth_',num2str(binwidth),'.mat'])
end


function [Pcon_inh,Pcon_exc,CV_inh,CV_exc] = connection_prob_CV_averaged(rin,a,trial_num_range,binwidth)

Pcon_inh_trial = nan(length(trial_num_range),1);
Pcon_exc_trial = nan(length(trial_num_range),1);
CV_inh_trial = nan(length(trial_num_range),1);
CV_exc_trial = nan(length(trial_num_range),1);
for j = 1:length(trial_num_range)
    trial_num = trial_num_range(j);
    fname = ['/home/zhang.chi9/research/logscale//network_load_at_numerical_capacity/rin_',num2str(rin),'a_',num2str(a),'TrialNum_',num2str(trial_num),'.mat'];
    if isfile(fname)
        load(fname)
        [Pcon_inh_trial(j),Pcon_exc_trial(j),CV_inh_trial(j),CV_exc_trial(j)] = connection_prob_CV_fit(W,binwidth);
    end
end
Pcon_inh = nanmean(Pcon_inh_trial);
Pcon_exc = nanmean(Pcon_exc_trial);
CV_inh = nanmean(CV_inh_trial);
CV_exc = nanmean(CV_exc_trial);
end
