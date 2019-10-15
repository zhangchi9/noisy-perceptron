clear
close all
rin = 2^-6;
rout = 2^-6;
mean_log =  log(0.2);
var_log = 0.3;
F = lognrnd(mean_log,var_log,1,10^8);

pool=F(rin<2*F.*(1-F) & F < 0.5);
%pool = R(R<0.3 & R > 0.02);
N = 10000;
Ninh = N*0.2;
Nexc = N-Ninh;
f = [sort(datasample(pool,Ninh)),sort(datasample(pool,N-Ninh))]';
%f = [0.05*ones(Ninh/2,1);0.05*ones(Ninh/2,1);0.05*ones(Nexc/4,1);0.05*ones(Nexc/2,1);0.2*ones(Nexc/4,1);];
%f = [linspace(0.02,0.4,Ninh),linspace(0.02,0.4,N-Ninh)]';
%x = (0:0.001:0.5)';y = lognpdf(x,mean_log,var_log);figure; plot(x,y),axis square
data = datasample(pool,10^7);
figure
num_bars = 50; %// specify number of bars
[n, x] = hist(data,num_bars); %// use two-output version of hist to get values
n_normalized = n/numel(data)/(x(2)-x(1)); %// normalize to unit area
%bar(x, n_normalized, 1); %// plot histogram (with unit-width bars)
hold on
plot(x, n_normalized, 'k'); %// plot line, in red (or change color)
axis square
box on

fout = 0.2%[linspace(0.02,0.4,100)]';
Pcon = zeros(N,length(fout));
Jmean = zeros(N,length(fout));
for k = 1:length(fout)
[~,~,Pcontmp,~,Jmeantmp,~] = theoretical_solution_heter(40,0,rin,rout,f,fout(k),70,'heter',N);
Pcon(:,k) = Pcontmp;
Jmean(:,k) = Jmeantmp;
end

figure, plot(repmat(f(1:Ninh),1,length(fout)),Pcon(1:Ninh,:),'r'), hold on 
plot(repmat(f((Ninh+1):N),1,length(fout)),Pcon((Ninh+1):N,:),'b')
title('Pcon1')
legend('Inh','Exc')
axis square

figure, plot(repmat(f(1:Ninh),1,length(fout)),Jmean(1:Ninh,:),'r'), hold on 
plot(repmat(f((Ninh+1):N),1,length(fout)),Jmean((Ninh+1):N,:),'b')
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