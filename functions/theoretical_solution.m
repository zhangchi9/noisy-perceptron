% This function will calculate the robust model theoretical_solution,
% f: presynaptic neuron firing rate, p presynaptic neuron error rate
% fout: postsynaptic neuron firing rate, pout presynaptic neuron error rate
% w: L1 norm, finh: inhibitory neuron percentage
% All the probability are independent of j

function [capacity,exitflag,Pconinh,Pconexc,CVinh,CVexc] = theoretical_solution(rou,beta,rin,f)

N = 1000;
%f = 0.2;
fout = 0.2;
w = 70;
finh = 0.2;
Ninh = N*finh;

rout = rin;
pout1 = rout/2./(1-fout);
pout2 = (1-fout)./fout.*pout1;
Cout1 = sqrt(2)*erfinv(1-2*pout1);
Cout2 = sqrt(2)*erfinv(1-2*pout2);

E = @(x) (1+erf(x))/2;
F = @(x) exp(-x.^2)./pi^0.5+x.*(1+erf(x));
D = @(x) x.*F(x)+E(x);

g=zeros(N,1);
g(1:Ninh) = -1;
g(Ninh+1:end) = 1;

p1 = rin/2./(1-f);
p2 = (1-f)./f.*p1;
Cj = (1-f).*p1.*(1-p1)+ f.*p2.*(1-p2);
Dj = f.*(1-f).*(1-p1-p2).^2;

options = optimset('Display','off','MaxIter',10^5,'MaxFunEvals',10000,'TolX',10^-8,'TolFun',10^-8);

S = @(x) [(1-fout)*F(x(1))-fout*F(x(2));...
    
((1-finh)*F(x(4))+finh*F(x(3)))*(rou^2 + w*beta*(f+Cj))-(Dj*4*(x(1)+x(2))/(Cout1 + Cout2)^2*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))) + Cj)*w^2*((x(3)-x(4))/w/f-(x(3)+x(4)));...

((1-finh)*F(x(4))-finh*F(x(3)))*(rou^2 + w*beta*(f+Cj))-(Dj*4*(x(1)+x(2))/(Cout1 + Cout2)^2*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))) + Cj)*w*((x(3)-x(4))/w/f-(x(3)+x(4)))/f;...

((1-finh)*D(x(4))+finh*D(x(3)))*(rou^2 + w*beta*(f+Cj))*(2*Dj*(x(1)+x(2))^2/(Cout1 + Cout2)^2 - Cj)-(Dj*4*(x(1)+x(2))/(Cout1 + Cout2)^2*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))) + Cj)^2*w^2*((x(3)-x(4))/w/f-(x(3)+x(4)))^2/2];


for j = 1: 5000
    j
    x0 = [0.1 1.5 0.5 -0.5]+0.1.*rand(1,4);
    [x,~,exitflag] = fsolve(S, x0, options);
    if exitflag == 1 && ((x(3)-x(4))/w/f-(x(3)+x(4)))*(x(1) + x(2))>0 && (x(1) + x(2))/(Cout1 + Cout2)>=0
        break;
    end
end
if  exitflag == 1
    capacity = 16*(x(1)+x(2))^2/(Cout1 + Cout2)^4*(fout*D(x(2))+(1-fout)*D(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)))^2*((1-finh)*D(x(4)) + finh*D(x(3)))*Dj^2/(Dj*4*(x(1)+x(2))/(Cout1 + Cout2)^2*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))) + Cj)^2;
    Pconexc = E(x(4));
    Pconinh = E(x(3));
    CVinh = sqrt(2*D(x(3))*E(x(3))/F(x(3))^2 - 1);
    CVexc = sqrt(2*D(x(4))*E(x(4))/F(x(4))^2 - 1);
end

