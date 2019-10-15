function [mean_inh,mean_exc,Pcon_inh,Pcon_exc,CV_inh,CV_exc,std_exc,std_inh] = connection_prob_CV_add_rectangle(W,thinh,thexc)
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