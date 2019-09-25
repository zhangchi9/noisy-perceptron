function structure_analysis_map_new_fit()
%close all
rin_ind_range = [-2:-0.5:-12];
a_range = [0.1,5:5:100];
% rin_ind_range = [-2.5];
% a_range = [50];
%trial_range = [1:100];
Pcon_inh = nan(length(rin_ind_range),length(a_range));
Pcon_exc = nan(length(rin_ind_range),length(a_range));
CV_inh = nan(length(rin_ind_range),length(a_range));
CV_exc = nan(length(rin_ind_range),length(a_range));

% RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
% scratch_dir = ['/home/zhang.chi9/matlabtmp/', num2str(randi(10^9))];
% mkdir(scratch_dir);
% pc = parcluster('local');
% pc.JobStorageLocation = scratch_dir;
% parpool(pc,pc.NumWorkers);
% par
for i = 1:length(rin_ind_range)
    tmp1 = nan(1,length(a_range));
    tmp2 = nan(1,length(a_range));
    tmp3 = nan(1,length(a_range));
    tmp4 = nan(1,length(a_range));
    for j = 1:length(a_range)
        [rin_ind_range(i),a_range(j)]
        %[tmp1(j),tmp2(j),tmp3(j),tmp4(j),fitobject_inh,fitobject_exc] = ...
        %connection_prob_CV_averaged(rin_ind_range(i),a_range(j));
        %thinh = 0 + 100*(j-1)/20;
        thinh = 50;
%         if i==1
%             thinh = 140;
%         end
        thexc = thinh;
        [tmp1(j),tmp2(j),tmp3(j),tmp4(j)] = connection_prob_CV_averaged(rin_ind_range(i),a_range(j),thinh,thexc);
    end
    Pcon_inh(i,:) = tmp1;
    Pcon_exc(i,:) = tmp2;
    CV_inh(i,:) = tmp3;
    CV_exc(i,:) = tmp4;
end
save(['fmincon_structure_map_load_at_numerical_capacity_th',num2str(thinh),'and_140.mat'])
end

function [Pcon_inh,Pcon_exc,CV_inh,CV_exc] = connection_prob_CV_averaged(rin,a,thinh,thexc)%[Pconinh,Pconexc,CVinh,CVexc,fitobject_inh,fitobject_exc] = connection_prob_CV_averaged(rin,a)

W_all = [];

for trial_ind = 1:100
    %filename = ['~/research/logscale/network_different_size/rin_',num2str(rin),'a_',num2str(a),'TrialNum_',num2str(trial_ind),'N_',num2str(400),'.mat'];
    filename = ['/home/zhang.chi9/research/logscale//network_load_at_numerical_capacity/rin_',num2str(rin),'a_',num2str(a),'TrialNum_',num2str(trial_ind),'.mat'];
    if exist(filename)
        load(filename)
        W_all = [W_all,W];
    end
end
%figure,hist(W_all(:),1000),xlim([-500,500]),ylim([0,20000]),title('combined weights')
% figure,hist(W_nonf(:),1000),xlim([-500,500]),ylim([0,20000]),title('non feasible weights')
% figure,hist(W_fea(:),1000),xlim([-500,500]),ylim([0,20000]),title('feasible weights')

%[c_exc,c_inh,s] = theoretical_solution_fit(a,0,rj,f);
% [Pconinh,CVinh,fitobject_inh] = Pcon_CV_fit(WI,true,c_init_inh,s_init_inh);
% [Pconexc,CVexc,fitobject_exc] = Pcon_CV_fit(WE,true,c_init_exc,s_init_exc);
%thexc = 15; thinh = thexc;
[~,~,Pcon_inh,Pcon_exc,CV_inh,CV_exc,~,~] = connection_prob_CV(W_all,thinh,thexc);


% [~,~,Pcon_inh,Pcon_exc,CV_inh,CV_exc,~,~] = connection_prob_CV(W_fea,thinh,thexc);
% 
% 
% [~,~,Pcon_inh,Pcon_exc,CV_inh,CV_exc,~,~] = connection_prob_CV(W_nonf,thinh,thexc);

end

function [mean_inh,mean_exc,Pcon_inh,Pcon_exc,CV_inh,CV_exc,std_exc,std_inh] = connection_prob_CV(W,thinh,thexc)
N=size(W,1);
WE = W(N*0.2+1 : end,: );
WI = W(1:N*0.2 , : );

WE_copy = WE(:);
WI_copy = abs(WI(:));

edges_I = 0:2:max(WI_copy);
bin_centers_I = (edges_I(1:end-1) + edges_I(2:end))/2;
N_counts_I = histcounts(WI_copy,edges_I);
%figure,plot(bin_centers_I,N_counts_I),title('inh')

edges_E = 0:2:max(WE_copy);
bin_centers_E = (edges_E(1:end-1) + edges_E(2:end))/2;
N_counts_E = histcounts(WE_copy,edges_E);
%figure,plot(bin_centers_E,N_counts_E),title('exc')

WI_copy = WI_copy(WI_copy>thinh);
WE_copy = WE_copy(WE_copy>thexc);

N_th_I = find(bin_centers_I<thinh,1,'last');
WI_copy = [WI_copy;thinh*rand(N_counts_I(N_th_I)*N_th_I,1)];
%figure,hist(WI_copy,bin_centers_I)

N_th_E = find(bin_centers_E<thexc,1,'last');
WE_copy = [WE_copy;thexc*rand(N_counts_E(N_th_E)*N_th_E,1)];
%figure,hist(WE_copy,bin_centers_E)

mean_exc = mean(WE_copy);
mean_inh = mean(WI_copy);

%Pcon_exc = sum(WE(:)>thexc)/numel(WE);
%Pcon_inh = sum(WI(:)<-thinh)/numel(WI);
Pcon_inh = numel(WI_copy(:))/numel(WI(:));
Pcon_exc = numel(WE_copy(:))/numel(WE(:));

%CV_exc = std(WE(WE(:)>thexc))/mean(WE(WE(:)>thexc));
%CV_inh = abs(std(WI(WI(:)<-thinh))/mean(WI(WI(:)<-thinh)));
CV_inh = abs(std(WI_copy(:)))/mean(WI_copy(:));
CV_exc = abs(std(WE_copy(:)))/mean(WE_copy(:));

std_exc = std(WE_copy);
std_inh = abs(std(WI_copy));
end

function [Pcon,CV,fitobject] = Pcon_CV_fit(W1,if_plot,c_init,s_init)
E = @(x) (1+erf(x))/2;
F = @(x) exp(-x.^2)./pi^0.5+x.*(1+erf(x));

J_av = s_init*F(-c_init/s_init/2^0.5)/sqrt(2)/E(-c_init/s_init/2^0.5);

n = numel(W1);
W1=abs(W1(:));
binwidth = 5;

% W1(W1<thr)=[];
bins = binwidth/2:binwidth:max(W1(:));
[y,~]=hist(W1(:),bins);
yden = y/sum(y)/binwidth;
lower_bound = [0,-inf,0,0];
upper_bound = [inf,500,500,10];
start_point = [0.1,0.01,100,1];
fit_options = fitoptions('Method','NonlinearLeastSquares','MaxIter',100000,'TolFun',10^-18,'TolX',10^-18,...
    'Lower',lower_bound,...
    'Upper',upper_bound,...
    'StartPoint',start_point);

%fitType = fittype('2/sqrt(2*pi)/s/(1-erf((thr+c)/sqrt(2)/s))*exp(-(x+c).^2/2/s^2)','problem','thr','options',fit_options);
%fitType = fittype('1/sqrt(2*pi)/s*exp(-(x/sqrt(2)/s+c).^2) + (1+erf(c))/sqrt(2*pi)/s0*exp(-x.^2/2/s0^2)','options',fit_options);
%fitType = fittype('1/sqrt(2*pi)/s*exp(-(x/sqrt(2)/s+c).^2) + A*exp(-x.^2/2/s0^2)','options',fit_options);
fitType = fittype('2/sqrt(2*pi)/s/(1-erf((thr+c)/sqrt(2)/s))*exp(-(x+c).^2/2/s^2)','options',fit_options);

[fitobject,gof,output] = fit(bins',yden',fitType);

% c = -10;
% s = 10;
% s0 = 1;
% fun_test = @(x) 1/sqrt(2*pi)/s*exp(-(x+c).^2/2/s^2) + (1+erf((c/s/sqrt(2))))/sqrt(2*pi)/s0*exp(-x.^2/2/s0^2);

new_bins=0.1:0.2:max(W1(:));
if if_plot
    figure,plot(bins,yden,'b-'),hold on
    plot(-1000:0.2:max(W1(:)),fitobject(-1000:0.2:max(W1(:))),'r-')
end

pp = @(x) 2/sqrt(2*pi)/fitobject.s/(1-erf((fitobject.c)/sqrt(2)/fitobject.s))*exp(-(x+fitobject.c).^2/2/fitobject.s^2);

% Pcon= sum(fitobject(new_bins))*(new_bins(2)-new_bins(1))*sum(y)/n;
% X=sum(fitobject(new_bins)'.*new_bins)/sum(fitobject(new_bins));
% X2=sum(fitobject(new_bins)'.*new_bins.^2)/sum(fitobject(new_bins));

% calculation based on fit
Pcon= E(-fitobject.c/fitobject.s/2^0.5);
X=sum(pp(new_bins).*new_bins)/sum(pp(new_bins));
X2=sum(pp(new_bins).*new_bins.^2)/sum(pp(new_bins));
CV=(X2-X^2)^0.5/X;

end

% 
% function [Pcon,CV,fitobject] = Pcon_CV_fit(W1,if_plot,c_init,s_init)
% E = @(x) (1+erf(x))/2;
% F = @(x) exp(-x.^2)./pi^0.5+x.*(1+erf(x));
% 
% J_av = s_init*F(-c_init/s_init/2^0.5)/sqrt(2)/E(-c_init/s_init/2^0.5);
% 
% n = numel(W1);
% W1=abs(W1(:));
% binwidth = J_av/20;
% thr = J_av/3+0;
% W1(W1<thr)=[];
% bins = thr+binwidth/2:binwidth:max(W1(:));
% [y,~]=hist(W1(:),bins);
% yden = y/sum(y)/binwidth;
% lower_bound = [-inf,0];
% upper_bound = [inf,s_init]*5;
% start_point = [c_init*0,s_init];
% fit_options = fitoptions('Method','NonlinearLeastSquares','MaxIter',100000,'TolFun',10^-12,'TolX',10^-12,...
%     'Lower',lower_bound,...
%     'Upper',upper_bound,...
%     'StartPoint',start_point);
% 
% fitType = fittype('2/sqrt(2*pi)/s/(1-erf((thr+c)/sqrt(2)/s))*exp(-(x+c).^2/2/s^2)','problem','thr','options',fit_options);
% [fitobject,gof,output] = fit(bins',yden',fitType,'problem',thr);
% R2 = gof.adjrsquare;
% 
% new_bins=0.1:0.2:max(W1(:));
% if if_plot
%     figure,plot(bins,yden,'b-'),hold on
%     plot(-3000:0.2:max(W1(:)),fitobject(-3000:0.2:max(W1(:))),'r-')
% end
% 
% pp = @(x) 2/sqrt(2*pi)/fitobject.s/(1-erf((fitobject.c)/sqrt(2)/fitobject.s))*exp(-(x+fitobject.c).^2/2/fitobject.s^2);
% 
% % Pcon= sum(fitobject(new_bins))*(new_bins(2)-new_bins(1))*sum(y)/n;
% % X=sum(fitobject(new_bins)'.*new_bins)/sum(fitobject(new_bins));
% % X2=sum(fitobject(new_bins)'.*new_bins.^2)/sum(fitobject(new_bins));
% 
% % calculation based on fit
% Pcon= E(-fitobject.c/fitobject.s/2^0.5);
% X=sum(pp(new_bins).*new_bins)/sum(pp(new_bins));
% X2=sum(pp(new_bins).*new_bins.^2)/sum(pp(new_bins));
% CV=(X2-X^2)^0.5/X;
% 
% end