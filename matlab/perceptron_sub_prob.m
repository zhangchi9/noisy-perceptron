function perceptron_sub_prob(N,alpha)

addpath('/home/zhang.chi9/research/logscale/code/functions')
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

a = 40;
rj = 2^-6;
rout = rj;
f = 0.2;
n_trials = 200;

m = round(alpha*N*0.1872);

W = nan(N,n_trials);
routmin = nan(1,n_trials);

create_parpool()
parfor i=1:n_trials
    i
    X = rand(N,m)<f; %Input patterns
    Xp = 2*(rand(1,m)<f)-1;
    while sum(Xp==-1) == m
        Xp = 2*(rand(1,m)<f)-1;
    end
    [W(:,i),routmin(i)] = perceptron_learning_with_noise(X,Xp,f,rj,a,rout);
end
delete(gcp)
save(sprintf('/scratch/zhang.chi9/perceptron/data/perceptron_suc_prob/N_%d_alpha_%f.mat',N,alpha))

