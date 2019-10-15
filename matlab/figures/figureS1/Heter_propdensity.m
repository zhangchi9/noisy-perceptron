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
[capacity,exitflag,Pcon,CV,Jmean,PropDens] = theoretical_solution_heter(40,0,rin,rout,f,'heter',N);

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

% figure,plot(Jmean(1:Ninh),Pcon(1:Ninh)), hold on,
% plot(Jmean((Ninh+1):N),Pcon((Ninh+1):N))
% xlabel('J')
% ylabel('Pcon')

% figure, plot(rin(1:Ninh),abs(CV(1:Ninh))), hold on 
% plot(rin((Ninh+1):N),CV((Ninh+1):N))
% title('CV1')
% legend('Inh','Exc')

% Jrange = (0:2000);
% 
% Ptot_exc = @(x) Propden_exc(x,PropDens,Pcon,N,Ninh);
% Ptot_inh = @(x) Propden_inh(x,PropDens,Pcon,Ninh);
% 
% figure, plot(Jrange,Ptot_exc(Jrange)),hold on, plot(-Jrange,Ptot_inh(-Jrange)),title('wdis1')
% axes('Position',[.65 .7 .2 .2])
% box on
% semilogy(-Jrange,Ptot_inh(-Jrange))
% hold on 
% semilogy(Jrange,Ptot_exc(Jrange))

%%
close all
rin = 2^-6;
rout = 2^-6;
mean_log =  log(0.2);
var_log = 1;
R = lognrnd(mean_log,var_log,1,50000);
pool = R(R<1);
N = 1000;
Ninh = N*0.2;
Nexc = N-Ninh;
f =[sort(datasample(pool,Ninh)),sort(datasample(pool,N-Ninh))]';
%f = [0.05*ones(Ninh/2,1);0.05*ones(Ninh/2,1);0.05*ones(Nexc/4,1);0.05*ones(Nexc/2,1);0.2*ones(Nexc/4,1);];
%f = [linspace(0.1,0.9,Ninh),linspace(0.1,0.9,N-Ninh)]';
x = (0:0.001:1)';y = lognpdf(x,mean_log,var_log);figure; plot(x,y),axis square
% f = [linspace(0.03,0.04,150),linspace(0.03,0.4,50),linspace(0.03,0.04,700),linspace(0.03,0.4,100)]';
% f = [sort(0.2+0.2.*randn(1,Ninh)), sort(0.2+0.2.*randn(1,800))]';
% f(f<0.03) = 0.03;
[capacity,exitflag,Pcon,CV,Jmean,PropDens] = theoretical_solution_heter(40,0,rin,rout,f,'heter',N);
%[~,~,Pcon_2,~,Jmean_2,~] = theoretical_solution_heter(40,0,rin,rout,[0.2*ones(N,1)],'heter',N);

figure, plot(f(1:Ninh),Pcon(1:Ninh),'-'), hold on 
plot(f((Ninh+1):N),Pcon((Ninh+1):N),'-')
title('Pcon2')
legend('Inh','Exc')
xlabel('f')
ylabel('Pcon')

figure, plot(f(1:Ninh),Jmean(1:Ninh),'-'), hold on 
plot(f((Ninh+1):N),Jmean((Ninh+1):N),'-')
title('non zeros connection weight average')
legend('Inh','Exc')
xlabel('f')
ylabel('weight')

% Jave = Jmean.*Pcon;
% figure, plot(f(1:Ninh),Jave(1:Ninh),'*-'), hold on 
% plot(f((Ninh+1):N),Jave((Ninh+1):N),'*-')
% Jave_2 = Jmean_2.*Pcon_2;
% plot([0,1],[Jave_2(1),Jave_2(1)])
% plot([0,1],[Jave_2(end),Jave_2(end)])
% title('J average')
% legend('Inh','Exc','Inh f 0.2','Exc f 0.2')
% xlabel('f')
% ylabel('J average')
% 
% current = Jmean.*Pcon.*f;
% figure, plot(f(1:Ninh),current(1:Ninh),'*-'), hold on 
% plot(f((Ninh+1):N),current((Ninh+1):N),'*-')
% current_2 = Jmean_2.*Pcon_2.*0.2;
% plot([0,1],[current_2(1),current_2(1)])
% plot([0,1],[current_2(end),current_2(end)])
% title('current')
% legend('Inh','Exc','Inh f 0.2','Exc f 0.2')
% xlabel('f')
% ylabel('current')
% 
% 
% figure,plot(Jmean(1:Ninh),Pcon(1:Ninh),'*-'), hold on,
% plot(Jmean((Ninh+1):N),Pcon((Ninh+1):N),'*-')
% xlabel('J')
% ylabel('Pcon')
% plot([0,0.8],[tmp(1),tmp(1)]*0.2)
% plot([0,0.8],[tmp(end),tmp(end)]*0.2)


% figure, plot(f(1:Ninh),abs(CV(1:Ninh))), hold on 
% plot(f((Ninh+1):N),CV((Ninh+1):N))
% title('CV2')
% legend('Inh','Exc')

% Jrange = (0:2000);
% 
% Ptot_exc = @(x) Propden_exc(x,PropDens,Pcon,N,Ninh);
% Ptot_inh = @(x) Propden_inh(x,PropDens,Pcon,Ninh);
% 
% figure, plot(Jrange,Ptot_exc(Jrange)),hold on, plot(-Jrange,Ptot_inh(-Jrange)),title('wdis2')
% axes('Position',[.65 .7 .2 .2])
% box on
% semilogy(Jrange,Ptot_exc(Jrange))
% hold on 
% semilogy(-Jrange,Ptot_inh(-Jrange))



