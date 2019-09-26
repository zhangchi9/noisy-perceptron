function Robust_heter_numerical_solution_rin_rout(rout,loadind)

N = 200;
Ninh = N * 0.2;
rj = [linspace(0.01,0.19,Ninh),linspace(0.01,0.19,N-Ninh)]';
rout = 10^-rout;

filename = (['N_',num2str(N),'_rout_',num2str(rout),'_load_',num2str(loadind)]);

%if ~(exist(filename,'file') == 2)
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

f = 0.2;
fout = f;
p1 = rj/2/(1-f);
p2 = (1-f)*p1/f;
memory_load = Robust_theoretical_solution_rin_rout(N,rj,rout);
m = round(N*memory_load*(1+0.1*(loadind-6)));
totol_training = 100;

W = nan(N,totol_training);
epsilonsum = nan(1,totol_training);
exitflag = nan(1,totol_training);

% scratch_dir = ['/home/zhang.chi9/matlabtmp/', num2str(randi(10^9))];
% mkdir(scratch_dir);
% pc = parcluster('local');
% pc.JobStorageLocation = scratch_dir;
% parpool(pc,56);
% par
for i = 1 : totol_training
    i
    X = rand(N,m)<f;
    Xp = 2*(rand(1,m)<fout)-1;
    X1 = (1 - p2*ones(1,m)).*X + (p1*ones(1,m)).*(1-X) ;
    X2 = (p2 * ones(1,m)).*(1 - p2 * ones(1,m)).*X + (p1 * ones(1,m)).*(1 ...
        - p1 * ones(1,m)).*(1-X);
    C = sqrt(2)*( erfinv(1-rout/fout)*(Xp+1)/2 + erfinv(1-rout/(1-fout))*...
        (1-(Xp+1)/2));
    [W,epsilonsum,exitflag] = find_W_single_trial(X1,X2,Xp,C);
    parsave([filename,'_trial_',num2str(i),'.mat'],X,Xp,rj,rout,memory_load,W,epsilonsum,exitflag)
end
delete(gcp)
end

function [W,epsilonsum,exitflag] = find_W_single_trial(X1,X2,Xp,C)
[N,m] = size(X1);
Ninh = N * 0.2;
w = 70;
g=[-ones(1,Ninh),ones(1,N-Ninh)];

delta = 10^-8;
opts = optimoptions(@fmincon,'Display','iter','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',false,'StepTolerance',delta,...
    'TolCon',delta,'TolX',delta,'TolFun',delta,'MaxIter',8000,...
    'MaxFunctionEvaluations',8000000);

fun = @(x) [zeros(1,N),ones(1,m)]*x;
AA = [-diag(g),zeros(N,m)];%sign constrainted
BB = [zeros(m,N),-eye(m,m)];
A = [AA;BB];
b = zeros(N+m,1);
Aeq = [g,zeros(1,m)];
beq = N*w;

nonlcon = @(x) nonlincon(x,X1,X2,Xp,C);
x0 = [rand(N,1)*70;100*ones(m,1)];

[SV,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[],[],nonlcon,opts);
W = SV(1:N);
epsilonsum = sum(SV(N+1:end));
end

function [y,yeq] = nonlincon(x,X1,X2,Xp,C)
[N,m] = size(X1);
J = x(1:N);
epsilon = x(N+1:end);
y = C./N.*sqrt(sum(( J.^2 * ones(1,m)).*X2)) - Xp.*( J'*X1/N - 1) - epsilon';
yeq = [];
end

function parsave(filename,X,Xp,rj,rout,memory_load,W,epsilonsum,exitflag)
save(filename)
end
