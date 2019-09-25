function Retrieval_prob_lenght_for_single_point()
rin = -6;
a = 40; 
prob = nan(1,20);
prob_std = nan(1,20);
retriveal_length = nan(1,20);
retriveal_length_std = nan(1,20);
for noise = 5:5:100
    index = round(noise/5);
    [prob(index),prob_std(index),retriveal_length(index),retriveal_length_std(index)] = retrieval(rin,a,noise);
end
figure,errorbar(5:5:100,prob,prob_std)
figure,errorbar(5:5:100,retriveal_length,retriveal_length_std)
save(['single_retriveal_length_prob,mat'])
end

function [prob,prob_std,retriveal_length,retriveal_length_std] = retrieval(rin,a,noise)
% rin = 0.02;
% a = 40; 
trial_num_range = 1:100;

retrieval_prob = nan(length(trial_num_range),100);
retrieval_length = nan(length(trial_num_range),100);
for j = 1:length(trial_num_range)
    j
    trial_num = trial_num_range(j);
    fname = ['/home/zhang.chi9/research/logscale/network_load_0.85/rin_',num2str(rin),'a_',num2str(a),'TrialNum_',num2str(trial_num),'loadfac_0.85.mat'];
    if isfile(fname)
        load(fname)
        if_retrieval = nan(1,100);
        len_retrieval = nan(1,100);
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
            tmpfind = find(hd>0.15,1,'first') - 1 ;
            if ~isempty(tmpfind)
                len_retrieval(ind)=tmpfind;
            else
                len_retrieval(ind) = m;
            end
        end
        %retrieval_prob(j) = nanmean(if_retrieval);
        retrieval_prob(j,:) = if_retrieval;
        retrieval_length(j,:) = len_retrieval;
    end
end
prob = nanmean(retrieval_prob(:));
prob_std = nanstd(retrieval_prob(:))/sqrt(nnz(~isnan(retrieval_prob)));
retriveal_length = nanmean(retrieval_length(:));
retriveal_length_std = nanstd(retrieval_length(:))/sqrt(nnz(~isnan(retrieval_length)));
end
