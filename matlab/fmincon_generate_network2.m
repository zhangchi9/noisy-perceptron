% this function is the save as the fmincon generate network, but it runs
% many trials in one code and load it at the numerical capacity

function fmincon_generate_network2(rj_ind,a,N)


rin_range = -2:-0.5:-12;
a_range = [0.1,5:5:100];

load('m_capacity.mat')

m = floor(m_capacity(find(rj_ind == rin_range),find(a == a_range))/400*N);

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

rj = 2^rj_ind;
%N = 200;
rout = rj;
f = 0.2;
fout = f;
p1 = rj/2/(1-f)*ones(N,1);
p2 = (1-f)*p1/f;

maxTrial = 100;
filename_check(rj_ind,a,maxTrial,N)

for TrialNum = 1 : maxTrial
    
    X =rand(N,m+1)<f; %Input patterns
    XX = X(:,1:end-1);
    
    W = nan(N,N);
    epsilonsum = nan(1,N);
    exitflag = nan(1,N);
    
    scratch_dir = ['/home/zhang.chi9/matlabtmp/', num2str(randi(10^9))];
    mkdir(scratch_dir);
    pc = parcluster('local');
    pc.JobStorageLocation = scratch_dir;
    parpool(pc,pc.NumWorkers);
    parfor i=1:N
        i
        Xp = 2*X(i,2:end)-1;
        X1 = (1 - p2*ones(1,m)).*XX + (p1*ones(1,m)).*(1-XX) ;
        X2 = (p2 * ones(1,m)).*(1 - p2 * ones(1,m)).*XX + (p1 * ones(1,m)).*(1 - p1 * ones(1,m)).*(1-XX);
        C = sqrt(2)*( erfinv(1-rout/fout)*(Xp+1)/2 + erfinv(1-rout/(1-fout))*(1-(Xp+1)/2));
        [W(:,i),epsilonsum(i),exitflag(i)] = find_W_single_trial(XX,X1,X2,Xp,C,a);
    end
    delete(gcp)
    save(filename_check(rj_ind,a,TrialNum,N))
end
end

function filename = filename_check(rj_ind,a,TrialNum,N)
filename = ['/home/zhang.chi9/research/logscale/network_different_size/rin_',num2str(rj_ind),'a_',num2str(a),'TrialNum_' num2str(TrialNum),'N_',num2str(N),'.mat'];
while isfile(filename)
    TrialNum = TrialNum + 1;
    if TrialNum > 100
        error('100 network has generated')
    end
    filename = ['/home/zhang.chi9/research/logscale/network_different_size/rin_',num2str(rj_ind),'a_',num2str(a),'TrialNum_' num2str(TrialNum),'N_',num2str(N),'.mat'];
end
end

function [W,epsilonsum,exitflag] = find_W_single_trial(X,X1,X2,Xp,C,a)

[N,m] = size(X1);
Ninh = N * 0.2;
w = 70;
g=[-ones(1,Ninh),ones(1,N-Ninh)];

delta = 10^-8;
opts = optimoptions(@fmincon,'Display','off','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',false,'StepTolerance',delta,...
    'TolCon',delta,'TolX',delta,'TolFun',delta,'MaxIter',3000,...
    'MaxFunctionEvaluations',3000000);


fun = @(x) [zeros(1,N),ones(1,m)]*x;
AA = [-diag(g),zeros(N,m)];%sign constrainted
BB = [zeros(m,N),-eye(m,m)];
A = [AA;BB];
b = zeros(N+m,1);
Aeq = [g,zeros(1,m)];
beq = N*w;


nonlcon = @(x) nonlincon(x,X,X1,X2,Xp,C,a);
x0 = [rand(N,1)*70;100*ones(m,1)];


[SV,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[],[],nonlcon,opts);
W = SV(1:N);
epsilonsum = sum(SV(N+1:end));

end

function [y,yeq] = nonlincon(x,X,X1,X2,Xp,C,a)
beta = 0;
[N,m] = size(X1);
Ninh = N * 0.2;
J = x(1:N);
epsilon = x(N+1:end);
g=[-ones(1,Ninh),ones(1,N-Ninh)];
y = C./N.*sqrt(N*a^2 + beta*(J.*g')'*X + sum(( J.^2 * ones(1,m)).*X2)) - Xp.*( J'*X1/N - 1) - epsilon';
yeq = [];

end

function parsave(filename,X,Xp,rj,rout,W,epsilonsum,exitflag)

save(filename)

end
