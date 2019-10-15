clear
load cv_Amplitude.mat
load cv_Amplitude_Feldmeyer.mat

x = [cv_Amplitude(:,1)]%;cv_Amplitude_Feldmeyer(:,1)];
y = [(cv_Amplitude(:,1).*cv_Amplitude(:,2)).^2]%;(cv_Amplitude_Feldmeyer(:,1).*cv_Amplitude_Feldmeyer(:,2)).^2];
[xData, yData] = prepareCurveData( x, y );
ind = find(y>0.8);
% Set up fittype and options.
ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0;
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'linear fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. x', 'lineaer fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel J
ylabel var
grid on

figure, plot(cv_Amplitude(:,1),cv_Amplitude(:,2),'.','MarkerSize',20)
hold on
plot(cv_Amplitude_Feldmeyer(:,1),cv_Amplitude_Feldmeyer(:,2),'.','MarkerSize',20)
xlabel J
ylabel CV

betaj = y./x;
[betaj,I] = sort(betaj);
J = x(I);
binsize = 0.1;
edges = 0:binsize:(max(betaj)+binsize);
edges_cen = (edges(1:end-1) + edges(2:end))/2;
Y = discretize(betaj,edges);

for cat = 1:length(edges_cen)
    
    catind = find(Y == cat);
    
    Jcat = J(catind);
    
    Jcatmean(cat) = mean(Jcat);
    
    Jcatvar(cat) = std(Jcat)/sqrt(length(Jcat));

end
figure, plot(betaj,J,'.','MarkerSize',20)
xlabel betaj
ylabel J
hold on 
errorbar(edges_cen,Jcatmean,Jcatvar)

% x = [cv_Amplitude(:,1);cv_Amplitude_Feldmeyer(:,1)];
% y = [(cv_Amplitude(:,1).*cv_Amplitude(:,2)).^2;(cv_Amplitude_Feldmeyer(:,1).*cv_Amplitude_Feldmeyer(:,2)).^2];
% 
% betaj = y./x;
% [betaj,I] = sort(betaj);
% J = x(I);
% binsize = 0.1;
% edges = 0:binsize:(max(betaj)+binsize);
% edges_cen = (edges(1:end-1) + edges(2:end))/2;
% Y = discretize(betaj,edges);
% 
% for cat = 1:length(edges_cen)
%     
%     catind = find(Y == cat);
%     
%     Jcat = J(catind);
%     
%     Jcatmean(cat) = mean(Jcat);
%     
%     Jcatvar(cat) = std(Jcat)/sqrt(length(Jcat));
% 
% end
% figure, plot(betaj,J,'.','MarkerSize',20)
% xlabel betaj
% ylabel J
% hold on 
% errorbar(edges_cen,Jcatmean,Jcatvar)
% 
% x = cv_Amplitude_Feldmeyer(:,1);
% y = (cv_Amplitude_Feldmeyer(:,1).*cv_Amplitude_Feldmeyer(:,2)).^2;
% 
% betaj = y./x;
% [betaj,I] = sort(betaj);
% J = x(I);
% binsize = 0.05;
% edges = 0:binsize:(max(betaj)+binsize);
% edges_cen = (edges(1:end-1) + edges(2:end))/2;
% Y = discretize(betaj,edges);
% Jcatmean = [];
% Jcatvar = [];
% for cat = 1:length(edges_cen)
%     
%     catind = find(Y == cat);
%     
%     Jcat = J(catind);
%     
%     Jcatmean(cat) = mean(Jcat);
%     
%     Jcatvar(cat) = std(Jcat)/sqrt(length(Jcat));
% 
% end
% figure, plot(betaj,J,'.','MarkerSize',20)
% xlabel betaj
% ylabel J
% hold on 
% errorbar(edges_cen,Jcatmean,Jcatvar)
% 
