% this function is the save as the fmincon generate network, but it runs
% many trials in one code and load it at the numerical capacity

function generate_network_N800(model,TrialNum,seperation)

N = 800;
rj_ind = -6;
beta_post = 40;

addpath('/home/zhang.chi9/research/logscale/code/functions')
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
rng(2019+TrialNum)


c = containers.Map;
% c('perceptron_20') = 1;
% c('fmincon_20') = 1;
c('perceptron_200') = 0.1243;
c('perceptron_400') = 0.1354;
c('perceptron_800') = 0.1466;
c('fmincon_200') = 0.1492;
c('fmincon_400') = 0.1639;
c('fmincon_800') = 0.1764;

m = round(c([model,'_',num2str(N)])*N);

rj = 2^rj_ind;
rout = rj;
f = 0.2;

X =rand(N,m+1)<f; %Input patterns

while any(mean(X,2) == 0)
    X =rand(N,m+1)<f;
end

XX = X(:,1:end-1);

W = nan(N,N);
epsilonsum = nan(1,N);
exitflag = nan(1,N);

create_parpool()
parfor i= 100*seperation +(1:100)
    i
    Xp = 2*X(i,2:end)-1;
    [W(:,i),epsilonsum(i),exitflag(i)] = numerical_solution(XX,Xp,f,rj,beta_post,rout,model);
end
delete(gcp)
save(['/scratch/zhang.chi9/perceptron/data/tmp_N_800/',model,'_',num2str(TrialNum),'_',num2str(seperation)])
end
