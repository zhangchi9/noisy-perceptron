function perceptron_sub_prob(N,alpha)

addpath('/home/zhang.chi9/research/logscale/code/functions')
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));


a = 40;
rj = 2^-6;
rout = rj;
f = 0.2;
fout = f;
n_trials = 200;

p1 = rj/2/(1-f)*ones(N,1);
p2 = (1-f)*p1/f;
m = round(alpha*N);
X =rand(N,m+1)<f; %Input patterns
XX = X(:,1:end-1);

W = nan(N,n_trials);
routmin = nan(1,N);

create_parpool
parfor i=1:n_trials
    i
    Xp = 2*X(i,2:end)-1;
    [W(:,i),routmin(i)] = perceptron_learning_with_noise(XX,Xp,f,rj,a,rout);
end
delete(gcp)
save(sprintf('/home/zhang.chi9/research/logscale/data/perceptron_suc_prob/N_%d_alpha_%f.mat',N,alpha))

