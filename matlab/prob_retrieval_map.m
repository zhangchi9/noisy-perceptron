function prob_retrieval_map()
addpath('/home/zhang.chi9/research/logscale/code/functions')
rin_range_ind =  -2:-0.5:-12;
a_range = [0.1,5:5:100];
trial_range = [1:100];

noise_added = [0.1,5:5:100];

retrieval_prob =  nan(length(rin_range_ind),length(a_range),length(noise_added));
retrieval_length =  nan(length(rin_range_ind),length(a_range),length(noise_added));
rout =  nan(length(rin_range_ind),length(a_range),length(noise_added));

create_parpool()
parfor i = 1:length(rin_range_ind)
    tmp = nan(length(a_range),length(noise_added));
    tmp2 = nan(length(a_range),length(noise_added));
    tmp3 = nan(length(a_range),length(noise_added));
    for j = 1:length(a_range)
        for k = 1:length(noise_added)
            [rin_range_ind(i),a_range(j),noise_added(k)]
            [tmp(j,k),tmp2(j,k),tmp3(j,k)] =  retrieval_intrinsic_and_pre(rin_range_ind(i),a_range(j),trial_range,noise_added(k));
        end
    end
    retrieval_prob(i,:,:) = tmp;
    retrieval_length(i,:,:) = tmp2;
    rout(i,:,:) = tmp3;
end
delete(gcp)
save(['/home/zhang.chi9/research/logscale/results/retrieval_prob_length_map_load_at_numerical_capacity_noise_all.mat'])
end

function [prob,length_retrieval,rout] = retrieval_intrinsic_and_pre(rin,a,trial_num_range,noise)
len_retrieval_ave = nan(length(trial_num_range),1);
retrieval_prob = nan(length(trial_num_range),1);
rout_j = nan(length(trial_num_range),1);
for j = 1:length(trial_num_range)
    trial_num = trial_num_range(j);
    fname = ['/scratch/zhang.chi9/perceptron/data/network_load_at_numerical_capacity/rin_',num2str(rin),'a_',num2str(a),'TrialNum_',num2str(trial_num),'.mat'];
    if isfile(fname)
        load(fname)
        nruns = 1000;
        if_retrieval = nan(1,nruns);
        rout_runs = nan(1,nruns);
        len_retrieval = nan(1,nruns);
        for ind = 1:nruns
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
            rout_runs(ind) = hd(end);
            if_retrieval(ind) = hd(end)< 0.15;
            tmpfind = find(hd>0.15,1,'first') - 1 ;
            if ~isempty(tmpfind)
                len_retrieval(ind)=tmpfind;
            else
                len_retrieval(ind) = m;
            end
        end
        rout_j(j) = nanmean(rout_runs); 
        retrieval_prob(j) = nanmean(if_retrieval);
        len_retrieval_ave(j) = mean(len_retrieval);
    end
end
rout = nanmean(rout_j);
prob = nanmean(retrieval_prob);
length_retrieval = nanmean(len_retrieval_ave);
end
