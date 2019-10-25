close all
clear
f = 0.2;
%rin = [linspace(0.01,0.19,100),linspace(0.19,0.2,100),linspace(0.01,0.19,400),linspace(0.19,0.2,400)]';
%rin = 0.25 - 2.^[(linspace(-12,-2,Ninh)),(linspace(-12,-2,800))]';
mean_log = log(2^-6);
var_log = 1;
R = lognrnd(mean_log,var_log,1,50000);
pool = R(R<0.25);
N = 10000;
Ninh = N*0.2;
rin =[sort(datasample(pool,Ninh)),sort(datasample(pool,N-Ninh))]';
% rin(1) = 0.25;
% rin(Ninh) = 0.01;
% rin(Ninh+1) = 0.25;
% rin(N) = 0.01;
x = (0:0.001:0.25)';y = lognpdf(x,mean_log,var_log);figure; plot(x,y),axis square,title('dis')
%figure,hist(rin,100)
rout = 2^-6;
[capacity,exitflag,Pcon,CV,Jmean,PropDens] = theoretical_solution_heter(40,0,rin,rout,f,N);

%rin = [(linspace(-12,-2,Ninh)),(linspace(-12,-2,800))]';
figure, plot(rin(1:Ninh),Pcon(1:Ninh)), hold on 
plot(rin((Ninh+1):N),Pcon((Ninh+1):N))
title('Pcon1')
legend('Inh','Exc')
axis square

figure, plot(rin(1:Ninh),Jmean(1:Ninh)), hold on 
plot(rin((Ninh+1):N),Jmean((Ninh+1):N))
title('Jmean1')
legend('Inh','Exc')
axis square