% this function will calculate the connection probablity and CV of
% conenction weiths by fitting the distribution of the connection weights
% distribution.

function [Pcon_inh,Pcon_exc,CV_inh,CV_exc] = connection_prob_CV_fit(W,binwidth)

N=size(W,1);
WE = W(N*0.2+1 : end,: );
WI = W(1:N*0.2 , : );

Pcon_exc = nan(10,10);
Pcon_inh = nan(10,10);
CV_exc = nan(10,10);
CV_inh = nan(10,10);
R2E = nan(10,10);
R2I = nan(10,10);

for i = 1:10
    for j = 1:10
        upper_bound = i*[1000,1000];
        start_point = j*[300,300];
        fit_options = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[-1000,0],...
            'Upper',upper_bound,...
            'StartPoint',start_point);
        try
            [Pcon_exc(i,j),CV_exc(i,j),R2E(i,j)] = Pcon_CV_fit(WE,fit_options,binwidth);
            [Pcon_inh(i,j),CV_inh(i,j),R2I(i,j)] = Pcon_CV_fit(WI,fit_options,binwidth);
        catch
            continue
        end
    end
end
[~,I] = max(R2E(:));
[a,b] = ind2sub(size(R2E),I);
Pcon_exc = Pcon_exc(a,b);
CV_exc = CV_exc(a,b);

[~,I] = max(R2I(:));
[a,b] = ind2sub(size(R2I),I);
Pcon_inh = Pcon_inh(a,b);
CV_inh = CV_inh(a,b);


end


function [Pcon,CV,R2] = Pcon_CV_fit(W1,fit_options,binwidth)
n = numel(W1);
W1=abs(W1(:));

thr=50;
W1(W1<thr)=[];


%binwidth = 20;

bins = thr+binwidth/2:binwidth:max(W1(:));

[y,~]=hist(W1(:),bins);

yden = y/sum(y)/binwidth;

% fit_options = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[-1000,0],...
%                'Upper',[1000,1000],...
%                'StartPoint',[400,300]);

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
