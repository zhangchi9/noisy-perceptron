% this function is the save as the fmincon generate network, but it runs
% many trials in one code and load it at the numerical capacity

function generate_network(N,rj_ind,beta_post,model,savepath)


rin_range = -2:-0.5:-12;
beta_post_range = [0.1,5:5:100];

load('m_capacity.mat')

m = floor(m_capacity(find(rj_ind == rin_range),find(beta_post == beta_post_range))/400*N);

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

rj = 2^rj_ind;
%N = 200;
rout = rj;
f = 0.2;
fout = f;
p1 = rj/2/(1-f)*ones(N,1);
p2 = (1-f)*p1/f;

maxTrial = 100;
filename_check(rj_ind,beta_post,maxTrial,N)

for TrialNum = 1 : maxTrial
    
    X =rand(N,m+1)<f; %Input patterns
    XX = X(:,1:end-1);
    
    W = nan(N,N);
    epsilonsum = nan(1,N);
    exitflag = nan(1,N);
    
    create_parpool()
    parfor i=1:N
        i
        Xp = 2*X(i,2:end)-1;
        [W(:,i),epsilonsum(i),exitflag(i)] = numerical_solution(XX,Xp,f,rj,beta_post,rout,model)
    end
    delete(gcp)
    save(filename_check(savepath,rj_ind,beta_post,TrialNum,N))
end
end