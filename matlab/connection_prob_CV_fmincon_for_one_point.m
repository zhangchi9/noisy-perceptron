
function connection_prob_CV_fmincon_for_one_point()

rin = -6;
a = 40;
trial_num_range = [1:100];
N = 400
Pcon_inh_trial = nan(length(trial_num_range),1);
Pcon_exc_trial = nan(length(trial_num_range),1);
CV_inh_trial = nan(length(trial_num_range),1);
CV_exc_trial = nan(length(trial_num_range),1);
for j = 1:length(trial_num_range)
    trial_num = trial_num_range(j);
    fname = ['/home/zhang.chi9/research/logscale/network_different_size/rin_',num2str(rin),'a_',num2str(a),'TrialNum_',num2str(trial_num),'N_',num2str(N),'.mat']
    if isfile(fname)
        load(fname)
        [Pcon_inh_trial(j),Pcon_exc_trial(j),CV_inh_trial(j),CV_exc_trial(j)] = connection_prob_CV_fit(W);
    end
end
Pcon_inh = nanmean(Pcon_inh_trial);
Pcon_exc = nanmean(Pcon_exc_trial);
CV_inh = nanmean(CV_inh_trial);
CV_exc = nanmean(CV_exc_trial);
%save('fmincon_structure_properties_for_one_point.mat')
end
