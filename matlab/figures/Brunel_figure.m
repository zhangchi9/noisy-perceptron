% This function will calculate the robust model theoretical_solution,
% f: presynaptic neuron firing rate, p presynaptic neuron error rate
% fout: postsynaptic neuron firing rate, pout presynaptic neuron error rate
% w: L1 norm, finh: inhibitory neuron percentage
% All the probability are independent of j

function [capacity,Pconexc] = Brunel_figure(rou,rout)

N = 1000;
f = 0.042;
fout = 0.37;
%w = 10;
finh = 0;
Ninh = N*finh;
beta = 0;
rin = 0;

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

S = @(x) [ (1-fout)*F(x(1))-fout*F(x(2));...

1/N*sum(f.*g.*sqrt(Dj)./(Cj + Dj*4*(x(1)+x(2))./(Cout1 + Cout2).^2*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)))).*F(-(2*f.*g*x(3) + beta*(f+Cj)*x(4))/2./sqrt(Dj))) - 2*x(4);...

1/N*sum(Dj./(Cj + Dj*4*(x(1)+x(2))./(Cout1 + Cout2).^2*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)))).^2.*(Cj/2-Dj*(x(1)+x(2)).^2/(Cout1 + Cout2).^2).*D(-(2*f.*g*x(3) + beta*(f+Cj)*x(4))/2./sqrt(Dj))) - rou.^2*x(4).^2 + 2*x(4)*x(3);...

-4*x(3)+4*rou.^2*x(4);];

for j = 1:5000
    j
    x0 = [0 1 0 0]+0.1.*rand(1,4);
    [x,~,exitflag] = fsolve(S, x0, options);
    if exitflag ==1 && x(4)*(x(1) + x(2)) > 0 && (x(1) + x(2))/(Cout1 + Cout2) >=0 && all(x(4)*(Cj + Dj*4*(x(1)+x(2))./(Cout1 + Cout2).^2*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)))) > 0)
        break;
    end
end

if  exitflag >0
    fac = 16*(x(1)+x(2)).^2/(Cout1+Cout2).^4*(fout*D(x(2))+(1-fout)*D(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))).^2;
    capacity = fac* 1/N*sum(Dj.^2./(Cj + Dj*4*(x(1)+x(2))/(Cout1 + Cout2).^2*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)))).^2.*D(-(2*f.*g*x(3) + beta*(f+Cj)*x(4))/2./sqrt(Dj)));
    Pcon = E(-(2*f.*g*x(3) + beta*(f+Cj)*x(4))/2./sqrt(Dj));
    Pconexc = Pcon(end);
% %     Jmean = sqrt(Dj)./(Cj + Dj*4*(x(1)+x(2))./(Cout1 + Cout2).^2*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))))/x(4)/2.*F(-(2*f.*g*x(3) + beta*(f+Cj)*x(4))/2./sqrt(Dj))./E(-(2*f.*g*x(3) + beta*(f+Cj)*x(4))/2./sqrt(Dj));
% %     CV = sqrt(2*D(-(2*f.*g*x(3) + beta*(f+Cj)*x(4))/2./sqrt(Dj)).*E(-(2*f.*g*x(3) + beta*(f+Cj)*x(4))/2./sqrt(Dj))./F(-(2*f.*g*x(3) + beta*(f+Cj)*x(4))/2./sqrt(Dj)).^2 - 1);
% %     % b = sqrt(2)*w*thetaj
% %     b = sqrt(Dj)./(Cj + Dj*4*(x(1)+x(2))./(Cout1 + Cout2).^2*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))))/x(4);
% %     c = -g.*(2*f.*g*x(3) + beta*(f+Cj)*x(4))/2./sqrt(Dj);
% %     % J is Jrange,b is sumj, c is pcon, d is Wag
% %     PropDensinh = @(J,i) exp(-(J./b(i) - c(i)).^2)./sqrt(pi)./b(i)./Pcon(i);
    
else
    capacity = NaN;
    Pconinh = NaN;
    Pconexc = NaN;
end


end