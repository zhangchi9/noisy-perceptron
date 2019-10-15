function [W,epsilonsum,exitflag] = fmincon_solution(X,Xp,f,rj,beta_post,rout)

p1 = rj/2/(1-f);
p2 = rj/2/f;
fout = f;

[N,m] = size(X);

X1 = (1 - p2*ones(1,m)).*X + (p1*ones(1,m)).*(1-X) ;
X2 = (p2 * ones(1,m)).*(1 - p2 * ones(1,m)).*X + (p1 * ones(1,m)).*(1 - p1 * ones(1,m)).*(1-X);
C = sqrt(2)*( erfinv(1-rout/fout)*(Xp+1)/2 + erfinv(1-rout/(1-fout))*(1-(Xp+1)/2));

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

nonlcon = @(x) nonlincon(x,X,X1,X2,Xp,C,beta_post);
x0 = [rand(N,1)*70;100*ones(m,1)];

[SV,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[],[],nonlcon,opts);
W = SV(1:N);
epsilonsum = sum(SV(N+1:end));

end

function [y,yeq] = nonlincon(x,X,X1,X2,Xp,C,beta_post)
beta = 0;
[N,m] = size(X1);
Ninh = N * 0.2;
J = x(1:N);
epsilon = x(N+1:end);
g=[-ones(1,Ninh),ones(1,N-Ninh)];
y = C./N.*sqrt(N*beta_post^2 + beta*(J.*g')'*X + sum(( J.^2 * ones(1,m)).*X2)) - Xp.*( J'*X1/N - 1) - epsilon';
yeq = [];
end