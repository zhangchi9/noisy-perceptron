function [Pcon_inh,Pcon_exc,CV_inh,CV_exc] = get_structure_different_size(file_obj,save_filename)

filename = {file_obj.name};
folder = {file_obj.folder};
W_all = [];
thinh=120;
thexc = 120;

for i = 1:length(filename)
    load([folder{i},'/',filename{i}])
    W_all = [W_all,W];
end

[~,~,Pcon_inh,Pcon_exc,CV_inh,CV_exc,~,~] = connection_prob_CV(W_all,thinh,thexc);

save(save_filename)
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