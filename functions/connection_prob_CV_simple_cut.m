function [mean_inh,mean_exc,Pcon_inh,Pcon_exc,CV_inh,CV_exc,std_exc,std_inh] = connection_prob_CV_simple_cut(W,thinh,thexc)
N=size(W,1);
WE = W(N*0.2+1 : end,: );
WI = W(1:N*0.2 , : );

mean_exc = mean(WE(WE(:)>thexc));
mean_inh = mean(WI(WI(:)<-thinh));

Pcon_exc = sum(WE(:)>thexc)/numel(WE);
Pcon_inh = sum(WI(:)<-thinh)/numel(WI);

CV_exc = std(WE(WE(:)>thexc))/mean(WE(WE(:)>thexc));
CV_inh = abs(std(WI(WI(:)<-thinh))/mean(WI(WI(:)<-thinh)));

std_exc = std(WE(WE(:)>thexc));
std_inh = abs(std(WI(WI(:)<-thinh)));

