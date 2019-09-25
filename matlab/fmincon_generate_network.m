function fmincon_generate_network(rj_ind,a,TrialNum,loadfac)

rj = 2^rj_ind;

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

N = 400;
rout = rj; 
f = 0.2;
fout = f;

alpha = theoretical_solution(a,0,rj,f,'homo');

p1 = rj/2/(1-f)*ones(N,1);
p2 = (1-f)*p1/f;
m = round(alpha*N*loadfac);


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
save(['/home/zhang.chi9/research/logscale/network_load_0.85/rin_',num2str(rj_ind),'a_',num2str(a),'TrialNum_' num2str(TrialNum) 'loadfac_' num2str(loadfac) '.mat'])

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
