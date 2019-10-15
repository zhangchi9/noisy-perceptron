function [prob,prob_std] = retrieval_intrinsic_and_pre(rin,a,trial_num_range,noise)

retrieval_prob = nan(length(trial_num_range),1);
for j = 1:length(trial_num_range)
    trial_num = trial_num_range(j);
    fname = ['rin_',num2str(rin),'a_',num2str(a),'TrialNum_',num2str(trial_num),'loadfac_0.85.mat'];
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
        end
        retrieval_prob(j) = nanmean(if_retrieval);
    end
end
prob = nanmean(retrieval_prob);
prob_std = nanstd(retrieval_prob);
end


