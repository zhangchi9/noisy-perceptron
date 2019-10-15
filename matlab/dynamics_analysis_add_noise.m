function dynamics_analysis_add_noise()
rin_range_ind = -2:-0.5:-12;
a_range = [0.1,5:5:100];
trial_range = 1:100;

CV_ISI = nan(length(rin_range_ind),length(a_range));
CV_ISI2 = nan(length(rin_range_ind),length(a_range));
SPKS_COR = nan(length(rin_range_ind),length(a_range));
COR_I = nan(length(rin_range_ind),length(a_range));

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
scratch_dir = ['/home/zhang.chi9/matlabtmp/', num2str(randi(10^9))];
mkdir(scratch_dir);
pc = parcluster('local');
pc.JobStorageLocation = scratch_dir;
parpool(pc,pc.NumWorkers);
parfor i = 1:length(rin_range_ind)
    tmp1 = nan(1,length(a_range));
    tmp2 = nan(1,length(a_range));
    tmp3 = nan(1,length(a_range));
    tmp4 = nan(1,length(a_range));
    for j = 1:length(a_range)
        [rin_range_ind(i),a_range(j)]
        [tmp1(j),tmp4(j),tmp2(j),tmp3(j)] = get_dynamics_add_noise(rin_range_ind(i),a_range(j),trial_range);
    end
    CV_ISI(i,:) = tmp1;
    CV_ISI2(i,:) = tmp4;
    SPKS_COR(i,:) = tmp2;
    COR_I(i,:) = tmp3;
end

save('dynamics_data_load_0.85_add_noise.mat')
end


