function [prob,length_retrieval_fraction,rout] = get_retrieval_prob_length(file_obj,noise_added,save_filename)

filename_list = {file_obj.name};
folder = {file_obj.folder};

len_retrieval_ave = nan(length(filename_list),1);
retrieval_prob = nan(length(filename_list),1);
rout_j = nan(length(filename_list),1);

for jj = 1:length(filename_list)
    load([folder{jj},'/',filename_list{jj}])
    nruns = 1000;
    if_retrieval = nan(1,nruns);
    rout_runs = nan(1,nruns);
    len_retrieval = nan(1,nruns);
    for ind = 1:nruns
        hd = zeros(1,size(X,2)-1);
        for i = 1 : (size(X,2)-1)
            if i == 1
                x = X(:,i);
            else
                x = y;
            end
            y = W'*x/N + randn(N,1)*noise_added/sqrt(N) > 1;
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
    rout_j(jj) = nanmean(rout_runs);
    retrieval_prob(jj) = nanmean(if_retrieval);
    len_retrieval_ave(jj) = mean(len_retrieval);
    
end
rout = nanmean(rout_j);
prob = nanmean(retrieval_prob);
length_retrieval_fraction = nanmean(len_retrieval_ave)/m;
if nargin == 3
save(save_filename)
end
end