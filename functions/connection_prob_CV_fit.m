% this function will calculate the connection probablity and CV of
% conenction weiths by fitting the distribution of the connection weights
% distribution. 

function [Pcon_inh,Pcon_exc,CV_inh,CV_exc] = connection_prob_CV_fit(W)

N=size(W,1);
WE = W(N*0.2+1 : end,: );
WI = W(1:N*0.2 , : );

[Pcon_exc,CV_exc,R2_exc] = Pcon_CV_fit(WE);
[Pcon_inh,CV_inh,R2_inh] = Pcon_CV_fit(WI);

end


function [Pcon,CV,R2] = Pcon_CV_fit(W1)
n = numel(W1);
W1=abs(W1(:));

thr=50;
W1(W1<thr)=[];


binwidth = 5;

bins = thr+binwidth/2:binwidth:max(W1(:));

[y,~]=hist(W1(:),bins);

yden = y/sum(y)/binwidth;

fit_options = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-1000,0],...
               'Upper',[1000,1000],...
               'StartPoint',[400,300]);
      
fitType = fittype('2/sqrt(2*pi)/s/(1-erf((thr+c)/sqrt(2)/s))*exp(-(x+c).^2/2/s^2)','problem','thr','options',fit_options);
[fitobject,gof,output] = fit(bins',yden',fitType,'problem',thr);
R2 = gof.adjrsquare;
% fit_options = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0,0],...
%                'Upper',[inf,inf,inf],...
%                'StartPoint',[500,500,100]);
% fitType = fittype('A*exp(-(x+c).^2/2/s^2)','options',fit_options);
% [fitobject,gof,output] = fit(bins',y',fitType)

%figure,plot(bins,yden,'b-'),hold on
new_bins=0.1:0.2:max(W1(:));
%plot(new_bins,fitobject(new_bins),'r-')
% Pcon=sum(fitobject(new_bins))*(new_bins(2)-new_bins(1))/(bins(2)-bins(1))/n;
% X=sum(fitobject(new_bins)'.*new_bins)/sum(fitobject(new_bins));
% X2=sum(fitobject(new_bins)'.*new_bins.^2)/sum(fitobject(new_bins));
% CV=(X2-X^2)^0.5/X;

Pcon= sum(fitobject(new_bins)*sum(y))*(new_bins(2)-new_bins(1))/n;
X=sum(fitobject(new_bins)'.*new_bins)/sum(fitobject(new_bins));
X2=sum(fitobject(new_bins)'.*new_bins.^2)/sum(fitobject(new_bins));
CV=(X2-X^2)^0.5/X;

end



% n = numel(W);
% thr=50;
% W=abs(W(:));
% W(W<thr)=[];
% bins = thr+2.5:5:max(W(:));
% 
% [y,~]=hist(W(:),bins);
% 
% fit_options = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0,-inf],...
%                'Upper',[inf,inf,inf],...
%                'StartPoint',[500,500,100]);
% fitType = fittype('b*exp(-(x+x0).^2/2/s2^2)','options',fit_options);
% [fitobject,gof,output] = fit(bins',y',fitType);
% R2 = gof.adjrsquare;
% %figure,plot(bins,y,'b-'),hold on
% new_bins=0.1:0.2:max(W(:));
% %plot(new_bins,fitobject(new_bins),'r-')
% Pcon=sum(fitobject(new_bins))*(new_bins(2)-new_bins(1))/(bins(2)-bins(1))/n;
% X=sum(fitobject(new_bins)'.*new_bins)/sum(fitobject(new_bins));
% X2=sum(fitobject(new_bins)'.*new_bins.^2)/sum(fitobject(new_bins));
% CV=(X2-X^2)^0.5/X;
% end
% 
