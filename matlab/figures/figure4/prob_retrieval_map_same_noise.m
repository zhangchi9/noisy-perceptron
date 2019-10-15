

function prob_retrieval_map_same_noise()
rin_range = [0.01:0.01:0.2];
a_range = [0.1,5:5:100];
trial_range = [1:50];

noise_added = [0.1,2,4,5:10:150];

retrieval_prob =  nan(length(rin_range),length(a_range),length(noise_added));
    
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
scratch_dir = ['/home/zhang.chi9/matlabtmp/', num2str(randi(10^9))];
mkdir(scratch_dir);
pc = parcluster('local');
pc.JobStorageLocation = scratch_dir;
parpool(pc,pc.NumWorkers);
parfor i = 1:length(rin_range)
    tmp = nan(length(a_range),length(noise_added));
    for j = 1:length(a_range)
        [rin_range(i),a_range(j)]
        for k = 1:length(noise_added)
            tmp(j,k) =  retrieval_intrinsic_and_pre(rin_range(i),a_range(j),trial_range,noise_added(k));
        end
    end
    retrieval_prob(i,:,:) = tmp;
end
delete(gcp)
save prob_retrieval_map_load85_same_noise.mat
end

function [prob,prob_std] = retrieval_intrinsic_and_pre(rin,a,trial_num_range,noise)

retrieval_prob = nan(length(trial_num_range),1);
for j = 1:length(trial_num_range)
    trial_num = trial_num_range(j);
    fname = ['/home/zhang.chi9/research/homo_network_load/network_load_0.85/rin_',num2str(rin),'a_',num2str(a),'TrialNum_',num2str(trial_num),'loadfac_0.85.mat'];
    if isfile(fname)
        load(fname)
        if_retrieval = nan(1,100);
        for ind = 1:100
            hd = zeros(1,size(X,2)-1);
            for i = 1 : (size(X,2)-1)
                if i == 1
                    x = X(:,i);
%                     p1 = rj/2/(1-f);
%                     p2 = rj/2/f;
%                     err1 = rand(N,1)<p1;    % 0 -> 1
%                     err2 = rand(N,1)<p2;    % 1 -> 0
%                     x(x==0) = x(x==0) + err1(x==0);
%                     x(x==1) = x(x==1) - err2(x==1);
%                     raster(:,i) = x;
                else
                    x = y;
                end
                y = W'*x/N + randn(N,1)*noise/sqrt(N) > 1;
                hd(i) = nnz(X(:,i+1)-y)/N;
            end
            if_retrieval(ind) = hd(end)< 0.15;
            final_rout = hd(end);
        end
        retrieval_prob(j) = nanmean(if_retrieval);
    end
end
prob = nanmean(retrieval_prob);
prob_std = nanstd(retrieval_prob);
end



