close all
clear
f = 0.2;
rin = [linspace(0.01,0.19,100),linspace(0.19,0.2,100),linspace(0.01,0.19,400),linspace(0.19,0.2,400)]';
rout = 0.002;
[capacity,exitflag,Pcon,CV,Jmean,PropDens] = theoretical_solution_heter(40,0,rin,rout,f,'heter');

figure, plot(rin(1:200),Pcon(1:200)), hold on 
plot(rin(201:1000),Pcon(201:1000))
title('Pcon1')
legend('Inh','Exc')

figure, plot(rin(1:200),Jmean(1:200)), hold on 
plot(rin(201:1000),Jmean(201:1000))
title('Jmean1')
legend('Inh','Exc')

figure, plot(rin(1:200),abs(CV(1:200))), hold on 
plot(rin(201:1000),CV(201:1000))
title('CV1')
legend('Inh','Exc')


Jrange = (0:2000);

Ptot_exc = @(x) Propden_exc(x,PropDens,Pcon);
Ptot_inh = @(x) Propden_inh(x,PropDens,Pcon);

figure, plot(Jrange,Ptot_exc(Jrange)),hold on, plot(-Jrange,Ptot_inh(-Jrange)),title('wdis1')
axes('Position',[.65 .7 .2 .2])
box on
semilogy(Jrange,Ptot_exc(Jrange))
hold on 
semilogy(-Jrange,Ptot_inh(-Jrange))
%%
%close all
rin = 0.02;
rout = 0.02;
f = [linspace(0.03,0.04,150),linspace(0.03,0.4,50),linspace(0.03,0.04,700),linspace(0.03,0.4,100)]';
% f = [sort(0.2+0.2.*randn(1,200)), sort(0.2+0.2.*randn(1,800))]';
% f(f<0.03) = 0.03;
[capacity,exitflag,Pcon,CV,Jmean,PropDens] = theoretical_solution_heter(40,0,rin,rout,f,'heter');

figure, plot(f(1:200),Pcon(1:200)), hold on 
plot(f(201:1000),Pcon(201:1000))
title('Pcon2')
legend('Inh','Exc')

figure, plot(f(1:200),Jmean(1:200)), hold on 
plot(f(201:1000),Jmean(201:1000))
title('Jmean2')
legend('Inh','Exc')

figure, plot(f(1:200),abs(CV(1:200))), hold on 
plot(f(201:1000),CV(201:1000))
title('CV2')
legend('Inh','Exc')

Jrange = (0:2000);

Ptot_exc = @(x) Propden_exc(x,PropDens,Pcon);
Ptot_inh = @(x) Propden_inh(x,PropDens,Pcon);

figure, plot(Jrange,Ptot_exc(Jrange)),hold on, plot(-Jrange,Ptot_inh(-Jrange)),title('wdis2')
axes('Position',[.65 .7 .2 .2])
box on
semilogy(Jrange,Ptot_exc(Jrange))
hold on 
semilogy(-Jrange,Ptot_inh(-Jrange))





