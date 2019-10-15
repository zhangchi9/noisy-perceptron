
function dynamics_analysis_add_noise()
rin_range = 0.01:0.01:0.2;
a_range = [0.1,5:5:100];
trial_range = 1:30;

CV_ISI = nan(length(rin_range),length(a_range));
CV_ISI2 = nan(length(rin_range),length(a_range));
SPKS_COR = nan(length(rin_range),length(a_range));
COR_I = nan(length(rin_range),length(a_range));

% RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
% scratch_dir = ['/home/zhang.chi9/matlabtmp/', num2str(randi(10^9))];
% mkdir(scratch_dir);
% pc = parcluster('local');
% pc.JobStorageLocation = scratch_dir;
% parpool(pc,pc.NumWorkers);
% par
for i = 1:length(rin_range)
    tmp1 = nan(1,length(a_range));
    tmp2 = nan(1,length(a_range));
    tmp3 = nan(1,length(a_range));
    tmp4 = nan(1,length(a_range));
    for j = 1:length(a_range)
        [rin_range(i),a_range(j)]
        [tmp1(j),tmp4(j),tmp2(j),tmp3(j)] = get_dynamics_add_noise(rin_range(i),a_range(j),trial_range);
    end
    CV_ISI(i,:) = tmp1;
    CV_ISI2(i,:) = tmp4;
    SPKS_COR(i,:) = tmp2;
    COR_I(i,:) = tmp3;
end

save('dynamics_data_load_0.85_add_noise.mat')
end

function [CV_ISI_mean,CV_ISI_mean2,SPKS_COR_mean,COR_I_mean] = get_dynamics_add_noise(rin,a,trial_num_range)

CV_ISI = nan(1,length(trial_num_range));
CV_ISI2 = nan(1,length(trial_num_range));
SPKS_COR = nan(1,length(trial_num_range));
COR_I = nan(1,length(trial_num_range));

for i = 1:length(trial_num_range)
    trial_num = trial_num_range(i);
    fname = ['/home/zhang.chi9/research/homo_network_load/network_load_0.85/rin_',num2str(rin),'a_',num2str(a),'TrialNum_',num2str(trial_num),'loadfac_0.85.mat'];
    if isfile(fname)
        load(fname)
        N = 400;
        n_run = 10;
        SPKS_COR_j = nan(1,n_run);
        COR_I_j = nan(1,n_run);
        CV_ISI_j = nan(1,n_run);
        CV_ISI_j2 = nan(1,n_run);
        for j = 1 : n_run
            [rin,a,i,j]
            [CV_ISI_all,CV_ISI_all2,SPKS_COR_j(j),~,~,~,COR_I_j(j)] = HammDist(W,a,N);
            CV_ISI_j(j) = nanmean(CV_ISI_all);
            CV_ISI_j2(j) = nanmean(CV_ISI_all2);
        end
        SPKS_COR(i) = nanmean(SPKS_COR_j);
        COR_I(i) = nanmean(COR_I_j);
        CV_ISI(i) = nanmean(CV_ISI_j);
        CV_ISI2(i) = nanmean(CV_ISI_j2);
    end
end
CV_ISI_mean = nanmean(CV_ISI);
CV_ISI_mean2 = nanmean(CV_ISI2);
SPKS_COR_mean = nanmean(SPKS_COR);
COR_I_mean = nanmean(COR_I);
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

