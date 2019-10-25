clear
%close all

R = 0.0156*ones(1,50000);

F = 0.2*ones(1,50000);

W = 70*ones(1,50000);

mean_log =  log(40);
var_log = 1;
BETA = lognrnd(mean_log,var_log,1,50000);
data = lognrnd(mean_log,var_log,1,10^7);
%BETA = 40*ones(1,50000);

ind=(R<2*F.*(1-F) & R>0 & F>0 & F<0.5);
R=R(ind);
F=F(ind);
W=W(ind);
BETA = BETA(ind);

N = 400;
Ninh = N*0.2;
rin = R(1:N)'; 
%rin = [linspace(0.01,0.25,Ninh),linspace(0.01,0.25,N-Ninh)]'; 
f = F(1:N)'; 
%f = [linspace(0.01,0.5,Ninh),linspace(0.01,0.5,N-Ninh)]'; 
w = W(1:N)';
beta = BETA(1:N)';

rout = rin;
fout = f;
Pcon = zeros(N,length(rout));
Jmean = zeros(N,length(rout));
for k = 1:length(rout)
[~,~,Pcontmp,~,Jmeantmp,~] = theoretical_solution_heter(beta(k),0,rin,rout(k),f,N);
Pcon(:,k) = Pcontmp;
Jmean(:,k) = Jmeantmp;
end


J_E_E = Jmean(Ninh+1:N,Ninh+1:N);
Pcon_E_E = Pcon(Ninh+1:N,Ninh+1:N);
J_I_I = Jmean(1:Ninh,1:Ninh);
Pcon_I_I = Pcon(1:Ninh,1:Ninh);
J_I_E = Jmean(1:Ninh,Ninh+1:N);
Pcon_I_E = Pcon(1:Ninh,Ninh+1:N);
J_E_I = Jmean(Ninh+1:N,1:Ninh);
Pcon_E_I = Pcon(Ninh+1:N,1:Ninh);

figure
num_bars = 5000; %// specify number of bars
[n, x] = hist(data,num_bars); %// use two-output version of hist to get values
n_normalized = n/numel(data)/(x(2)-x(1)); %// normalize to unit area
%bar(x, n_normalized, 1); %// plot histogram (with unit-width bars)
hold on
plot(x, n_normalized, 'k'); %// plot line, in red (or change color)
axis square
xlim([0,500])
box on

[beta_sorted,Ind] = sort(beta);
figure,hold on,plot(beta_sorted,Pcon(1,Ind),'r'),plot(beta_sorted,Pcon(end,Ind),'b')
box on 
axis square

figure,hold on,plot(beta_sorted,Jmean(1,Ind),'r'),plot(beta_sorted,Jmean(end,Ind),'b')
box on 
axis square
