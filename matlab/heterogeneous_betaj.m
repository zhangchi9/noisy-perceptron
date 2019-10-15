% This function will calculate the robust model theoretical_solution,
% f: presynaptic neuron firing rate, p presynaptic neuron error rate
% fout: postsynaptic neuron firing rate, pout presynaptic neuron error rate
% w: L1 norm, finh: inhibitory neuron percentage
% All the probability are independent of j

clear
N = 1000;
finh = 0.2;
%f = 0.2*rand(N,1)+0.1;
% f = normrnd(0.2,0.05,[N,1]);
% f(f<=0) = 0.2;
Ninh = N*finh;
finh = 0.2;
f = 0.2*ones(N,1);
betaj = [linspace(0.1,16,Ninh),linspace(0.1,16,N-Ninh)]';
%betaj = 0.1*normrnd(100,5,[N,1]);
%betaj(betaj<=0) = 10;
%betaj = exp(betaj);
fout = 0.2;
w = 70;
finh = 0.2;


error_out = 0.01;
pout1 = error_out/2./(1-fout);
pout2 = (1-fout)./fout.*pout1;
Cout1 = sqrt(2)*erfinv(1-2*pout1);
Cout2 = sqrt(2)*erfinv(1-2*pout2);

E = @(x) (1+erf(x))/2;
F = @(x) exp(-x.^2)./pi^0.5+x.*(1+erf(x));
D = @(x) x.*F(x)+E(x);

g=zeros(N,1);
g(1:Ninh) = -1;
g(Ninh+1:end) = 1;

options = optimset('Display','off','MaxIter',10^5,'MaxFunEvals',10000,'TolX',10^-8,'TolFun',10^-8);

Dj = f.*(1-f);

S = @(x) [fout*F(x(2))-(1-fout)*F(x(1)); ...
    sum(F(-(x(3) + 2*x(4)*f.*g + x(5)*f.*betaj)/2./sqrt(Dj))./sqrt(Dj))-8*N*w*x(5)*(x(1)+x(2))/(Cout1+Cout2)^2*...
    (fout*E(x(2))+(1-fout)*E(x(1)))/(fout*F(x(2))+(1-fout)*F(x(1))); ...
    sum(F(-(x(3) + 2*x(4)*f.*g + x(5)*f.*betaj)/2./sqrt(Dj)).*f.*g./sqrt(Dj))-8*N*x(5)*(x(1)+x(2))/(Cout1+Cout2)^2*...
    (fout*E(x(2))+(1-fout)*E(x(1)))/(fout*F(x(2))+(1-fout)*F(x(1)));...
    sum(F(-(x(3) + 2*x(4)*f.*g + x(5)*f.*betaj)/2./sqrt(Dj)).*f.*betaj./sqrt(Dj))-8*N*(w*x(3)+2*x(4))*(x(1)+x(2))/(Cout1+Cout2)^2*...
    (fout*E(x(2))+(1-fout)*E(x(1)))/(fout*F(x(2))+(1-fout)*F(x(1)));...
    sum(D(-(x(3) + 2*x(4)*f.*g + x(5)*f.*betaj)/2./sqrt(Dj)))-16*N*x(5)*(w*x(3)+2*x(4))/(Cout1+Cout2)^2*...
    (fout*E(x(2))+(1-fout)*E(x(1)))^2/(fout*F(x(2))+(1-fout)*F(x(1)))^2;];

j = 1;
while j < 500
    j
    x0 = [-1 1 0 1 0]+0.1.*rand(1,5);
    [x,~,exitflag] = fsolve(S, x0, options);
    if exitflag >0 && (x(1) + x(2))* x(5) > 0 &&(x(1) + x(2))/(Cout1 + Cout2) >=0
        break;
    end
    j = j + 1;
end
if  exitflag >0
     fac = 16/(Cout1 + Cout2)^2*(fout*D(x(2))+(1-fout)*D(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)))^2;
     capacity = fac* (w*x(3)+2*x(4))*x(5);
     Pcon = E(-(x(3) + 2*x(4)*f.*g + x(5)*f.*betaj)/2./sqrt(Dj));
     sumj = (Cout1+Cout2)^2/4/(x(1)+x(2))./sqrt(2*f.*(1-f))*(fout*F(x(2))+(1-fout)*F(x(1)))/(fout*E(x(2))+(1-fout)*E(x(1)));
     Jmean = g.*sumj/sqrt(2).*F(-(x(3) + 2*x(4)*f.*g + x(5)*f.*betaj)/2./sqrt(Dj))./E(-(x(3) + 2*x(4)*f.*g + x(5)*f.*betaj)/2./sqrt(Dj));
%     CV = g.*sqrt(2*D(-(x(3)+2*x(4)*f*g)/2/sqrt(f*(1-f))).*E(-(x(3)+2*x(4)*f*g)/2/sqrt(f*(1-f)))...
%         ./(F(-(x(3)+2*x(4)*f*g)/2/sqrt(f*(1-f)))).^2 - 1);
%     CVinh = -CV(1);
%     CVexc = CV(end);
end

capacity

figure
plot(betaj(1:Ninh), Pcon(1:Ninh),'.')
xlabel('betaj')
ylabel('Pcon')
title('Pcon inh')
figure
plot(betaj(Ninh+1:end), Pcon(Ninh+1:end),'.')
xlabel('betaj')
ylabel('Pcon')
title('Pcon exc')

figure
plot(betaj(1:Ninh), -Jmean(1:Ninh),'.')
xlabel('betaj')
ylabel('Jave')
title('Jave inh')
figure
plot(betaj(Ninh+1:end), Jmean(Ninh+1:end),'.')
xlabel('betaj')
ylabel('Jave')
title('Jave exc')


