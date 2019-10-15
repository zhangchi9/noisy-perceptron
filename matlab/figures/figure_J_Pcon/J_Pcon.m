close all
clear
f = 0.2;
fout = 0.2;
%rin = [linspace(0.01,0.19,100),linspace(0.19,0.2,100),linspace(0.01,0.19,400),linspace(0.19,0.2,400)]';
%rin = 0.25 - 2.^[(linspace(-12,-2,Ninh)),(linspace(-12,-2,800))]';
mean_log = log(2^-6);
var_log = 0.7;
%R = lognrnd(mean_log,var_log,1,50000);
R = randn(1,50000)*0.005*4 + 0.0156/2;
%R = rand(1,50000)*0.0156*2;
pool = R(R<0.25 & R > 0);
N = 800;
Ninh = N*0.2;
%rin =[sort(datasample(pool,Ninh)),sort(datasample(pool,N-Ninh))]';
rin = 2^-6;
% rin(1) = 0.25;
% rin(Ninh) = 0.01;
% rin(Ninh+1) = 0.25;
% rin(N) = 0.01;
x = (0:0.001:0.25)';y = lognpdf(x,mean_log,var_log);figure; plot(x,y),axis square,title('dis')
%figure,hist(rin,100)
rout = [sort(datasample(pool,100))]';
%rout = 2^-6;
Pcon = zeros(N,length(rout));
Jmean = zeros(N,length(rout));
for k = 1:length(rout)
[~,~,Pcontmp,~,Jmeantmp,~] = theoretical_solution_heter(40,0,rin,rout(k),f,fout,'heter',N);
Pcon(:,k) = Pcontmp;
Jmean(:,k) = Jmeantmp;
end


%rin = [(linspace(-12,-2,Ninh)),(linspace(-12,-2,800))]';
figure, plot(repmat(rin(1:Ninh),1,length(rout)),Pcon(1:Ninh,:),'.'), hold on 
plot(repmat(rin((Ninh+1):N),1,length(rout)),Pcon((Ninh+1):N,:),'.')
title('Pcon1')
legend('Inh','Exc')
axis square

figure, plot(repmat(rin(1:Ninh),1,length(rout)),Jmean(1:Ninh,:),'.'), hold on 
plot(repmat(rin((Ninh+1):N),1,length(rout)),Jmean((Ninh+1):N,:),'.')
title('Jmean1')
legend('Inh','Exc')
axis square

figure, plot(Jmean(1:Ninh,:),Pcon(1:Ninh,:),'.')
title('Inh')
xlabel('J')
ylabel('Pcon')
xlim([0,max(Jmean(1:Ninh,:))])

figure, plot(Jmean((Ninh+1):N,:), Pcon((Ninh+1):N,:),'.')
title('Exc')
xlabel('J')
ylabel('Pcon')
xlim([0,max(Jmean(Ninh+1:N,:))])
%%
rin = 2^-6;
rout = 2^-6;
mean_log =  log(0.1);
var_log = 0.5;
%R = lognrnd(mean_log,var_log,1,50000);
R = rand(1,50000)*0.2 + 0.1;
pool = R(R<0.3 & R > 0.02);
N = 1000;
Ninh = N*0.2;
Nexc = N-Ninh;
f =[sort(datasample(pool,Ninh)),sort(datasample(pool,N-Ninh))]';
%f = [0.05*ones(Ninh/2,1);0.05*ones(Ninh/2,1);0.05*ones(Nexc/4,1);0.05*ones(Nexc/2,1);0.2*ones(Nexc/4,1);];
%f = [linspace(0.02,0.4,Ninh),linspace(0.02,0.4,N-Ninh)]';
x = (0:0.001:0.5)';y = lognpdf(x,mean_log,var_log);figure; plot(x,y),axis square

fout = 0.2%[linspace(0.02,0.4,100)]';
Pcon = zeros(N,length(fout));
Jmean = zeros(N,length(fout));
for k = 1:length(fout)
[~,~,Pcontmp,~,Jmeantmp,~] = theoretical_solution_heter(40,0,rin,rout,f,fout(k),'heter',N);
Pcon(:,k) = Pcontmp;
Jmean(:,k) = Jmeantmp;
end

figure, plot(repmat(f(1:Ninh),1,length(fout)),Pcon(1:Ninh,:),'.'), hold on 
plot(repmat(f((Ninh+1):N),1,length(fout)),Pcon((Ninh+1):N,:),'.')
title('Pcon1')
legend('Inh','Exc')
axis square

figure, plot(repmat(f(1:Ninh),1,length(fout)),Jmean(1:Ninh,:),'.'), hold on 
plot(repmat(f((Ninh+1):N),1,length(fout)),Jmean((Ninh+1):N,:),'.')
title('Jmean1')
legend('Inh','Exc')
axis square

figure, plot(Jmean(1:Ninh,:),Pcon(1:Ninh,:),'.')
title('Inh')
xlabel('J')
ylabel('Pcon')

figure, plot(Jmean((Ninh+1):N,:),Pcon((Ninh+1):N,:),'.')
title('Exc')
xlabel('J')
ylabel('Pcon')

%%
w = 20:5:40;
Jexc = zeros(1,length(w));
Jinh = zeros(1,length(w));
Pconinh = zeros(1,length(w));
Pconexc = zeros(1,length(w));
for i = 1:length(w)
    [capacity,exitflag,Pconinh(i),Pconexc(i),CVinh,CVexc,Jexc(i),Jinh(i)] = theoretical_solution(w(i),0,2^-6,0.2,'homo');
end

figure, plot(w,Pconinh,'-'), hold on 
plot(w,Pconexc,'-')
title('Pcon3')
legend('Inh','Exc')
xlabel('beta')
ylabel('Pcon')

figure, plot(w,Jinh,'-'), hold on 
plot(w,Jexc,'-')
title('non zeros connection weight average')
legend('Inh','Exc')
xlabel('beta')
ylabel('weight')

figure, plot(Pconinh,Jinh)
title('Inh')
xlabel('Pcon')
ylabel('J')

figure, plot(Pconexc,Jexc)
title('Exc')
xlabel('Pcon')
ylabel('J')

%%
w = 30:5:300;
Jexc = zeros(1,length(w));
Jinh = zeros(1,length(w));
Pconinh = zeros(1,length(w));
Pconexc = zeros(1,length(w));
for i = 1:length(w)
    [capacity,exitflag,Pconinh(i),Pconexc(i),CVinh,CVexc,Jexc(i),Jinh(i)] = theoretical_solution(w(i),40,0,2^-6,0.2,'homo');
end

figure, plot(w,Pconinh,'-'), hold on 
plot(w,Pconexc,'-')
title('Pcon3')
legend('Inh','Exc')
xlabel('beta')
ylabel('Pcon')

figure, plot(w,Jinh,'-'), hold on 
plot(w,Jexc,'-')
title('non zeros connection weight average')
legend('Inh','Exc')
xlabel('beta')
ylabel('weight')

figure, plot(Pconinh,Jinh)
title('Inh')
xlabel('Pcon')
ylabel('J')

figure, plot(Pconexc,Jexc)
title('Exc')
xlabel('Pcon')
ylabel('J')