clear 
clc

addpath('/home/zhang.chi9/research/logscale/code/functions')
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

N = 200;
a = 40;
rj = 2^-6;
rout = rj;
f = 0.2;
fout = f;
alpha = 0.1;

p1 = rj/2/(1-f)*ones(N,1);
p2 = (1-f)*p1/f;
m = round(alpha*N);
X =rand(N,m+1)<f; %Input patterns
XX = X(:,1:end-1);

W = nan(N,N);
routmin = nan(1,N);

% create_parpool
% par
for i=1:N
    i
    Xp = 2*X(i,2:end)-1;
    [W(:,i),routmin(i)] = perceptron_learning_with_noise(XX,Xp,f,rj,a,rout);
end
delete(gcp)
save(sprintf('/home/zhang.chi9/research/logscale/perceptron_capacity/N_%d_TrialNum_%d.mat',N,TrialNum))

