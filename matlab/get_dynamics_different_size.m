function [CV_ISI,CV_ISI2,SPKS_COR,COR_I] = get_dynamics_different_size(file_obj,save_filename)

filename_list = {file_obj.name};
folder = {file_obj.folder};

CV_ISI = nan(1,length(filename_list));
CV_ISI2 = nan(1,length(filename_list));
SPKS_COR = nan(1,length(filename_list));
COR_I = nan(1,length(filename_list));

for i = 1:length(filename_list)
    load([folder{i},'/',filename_list{i}])
    
    n_run = 100;
    SPKS_COR_j = nan(1,n_run);
    COR_I_j = nan(1,n_run);
    CV_ISI_j = nan(1,n_run);
    CV_ISI_j2 = nan(1,n_run);
    for j = 1 : n_run
        [i,j]
        if exist('a','var')
            beta_post = a;
        end
        [CV_ISI_all,CV_ISI_all2,SPKS_COR_j(j),~,~,~,COR_I_j(j)] = HammDist(W,beta_post,N);
        CV_ISI_j(j) = nanmean(CV_ISI_all);
        CV_ISI_j2(j) = nanmean(CV_ISI_all2);
    end
    SPKS_COR(i) = nanmean(SPKS_COR_j);
    COR_I(i) = nanmean(COR_I_j);
    CV_ISI(i) = nanmean(CV_ISI_j);
    CV_ISI2(i) = nanmean(CV_ISI_j2);
    
end
% CV_ISI_mean = nanmean(CV_ISI);
% CV_ISI_mean2 = nanmean(CV_ISI2);
% SPKS_COR_mean = nanmean(SPKS_COR);
% COR_I_mean = nanmean(COR_I);

if nargin == 2
save(save_filename)
end
end


function [CV_ISI,CV_ISI2,SPKS_COR,IE,II,SDI,COR_I] = HammDist(W,a,N)
W = W';
CV_ISI = nan(N,1);
CV_ISI2 = nan(N,1);
SPKS_COR = nan;
IE = nan;
II = nan;
SDI = nan;
COR_I = nan;

Nruns = N;
S0 = rand(N,1) < 0.2;
S(:,1) = S0;
S2 (:,1) = S0;

Ninh = 0.2*N;
I_synE = nan(N,Nruns);
I_synI = nan(N,Nruns);

for i = 1:Nruns-1
    I_ext =0;%randn(N,1)*a/sqrt(N);
    S(:,i+1) = (W*S(:,i)./N -1 + I_ext)>0;
    
    S2(:,i+1) = (W*S(:,i)./N -1 + randn(N,1)*a/sqrt(N))>0;
    
    I_synE(:,i) = W(:,Ninh+1:end)*S(Ninh+1:end,i)/N;
    I_synI(:,i) = W(:,1:Ninh)*S(1:Ninh,i)/N;
    
end

idx = find(sum(S)>0);
if(length(idx)>5) %check if the spike raster has longer than 5 steps transient dynamics
    
    SDI = mean(std(I_synE(:,1:idx(end)-1) + I_synI(:,1:idx(end)-1),[],2)); %std across neurons
    IE =  mean(mean(I_synE(:,1:idx(end)-1)));
    II = mean(mean(I_synI(:,1:idx(end)-1)));
    
    tmpcorr = zeros(N,1);
    for n = 1:N
        tmp_corr = corrcoef(I_synE(n,1:idx(end)-1),-I_synI(n,1:idx(end)-1));
        tmpcorr(n) = tmp_corr(1,2);
    end
    COR_I = mean(tmpcorr); %correlation between excitatory and inhibitory currents
    
    %spike raster
    tmp_spikecorr = corrcoef(S(:,1:idx(end))');
    tmp_spikecorr(tmp_spikecorr==1)=nan;
    SPKS_COR = nanmean(tmp_spikecorr(:));
    
    for ll=1:N
        idx1 = find(S(ll,:)==1);
        ISI = diff(idx1);
        
        CV_ISI(ll) = std(ISI)/mean(ISI);
    end
end

tmpfind = find(mean(S)>0,1,'last')+1;
if ~isempty(tmpfind)
    S2(:,tmpfind:end) = 0;
end

idx = find(sum(S2)>0);
if(length(idx)>5) %check if the spike raster has longer than 5 steps transient dynamics
    
    SDI = mean(std(I_synE(:,1:idx(end)-1) + I_synI(:,1:idx(end)-1),[],2)); %std across neurons
    IE =  mean(mean(I_synE(:,1:idx(end)-1)));
    II = mean(mean(I_synI(:,1:idx(end)-1)));
    
    tmpcorr = zeros(N,1);
    for n = 1:N
        tmp_corr = corrcoef(I_synE(n,1:idx(end)-1),-I_synI(n,1:idx(end)-1));
        tmpcorr(n) = tmp_corr(1,2);
    end
    COR_I = mean(tmpcorr); %correlation between excitatory and inhibitory currents
    
    %spike raster
    tmp_spikecorr = corrcoef(S2(:,1:idx(end))');
    tmp_spikecorr(tmp_spikecorr==1)=nan;
    SPKS_COR = nanmean(tmp_spikecorr(:));
    
    for ll=1:N
        idx1 = find(S2(ll,:)==1);
        ISI = diff(idx1);
        
        CV_ISI2(ll) = std(ISI)/mean(ISI);
    end
end

end


