clear
f = 0.2;
fout = 0.2;

mean_log = log(2^-6);
var_log = 0.7;
R = lognrnd(mean_log,var_log,1,50000);
%R = randn(1,50000)*0.003 + 0.0156/2;
%R = rand(1,50000)*0.0156*2;
pool = R(R<2*f.*(1-f) & R>0);
N = 800;
Ninh = N*0.2;
rin =[sort(datasample(pool,Ninh)),sort(datasample(pool,N-Ninh))]';

% rin(1) = 0.25;
% rin(Ninh) = 0.01;
% rin(Ninh+1) = 0.25;
% rin(N) = 0.01;
figure,hist(rin,100)
%rout = [sort(datasample(pool,100))]';
rout = 2^-6;
Pcon = zeros(N,length(rout));
Jmean = zeros(N,length(rout));
for k = 1:length(rout)
[~,~,Pcontmp,~,Jmeantmp,~] = theoretical_solution_heter(40,0,rin,rout(k),f,fout,70,'heter',N);
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