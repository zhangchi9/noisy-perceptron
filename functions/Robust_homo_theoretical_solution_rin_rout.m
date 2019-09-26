% This function will calculate the robust model theoretical_solution,
% f: presynaptic neuron firing rate, p presynaptic neuron error rate
% fout: postsynaptic neuron firing rate, pout presynaptic neuron error rate
% w: L1 norm, finh: inhibitory neuron percentage
% All the probability are independent of j

function [capacity,exitflag,Pconinh,Pconexc,CVinh,CVexc,JmeanExc,JmeanInh,kappa] = Robust_homo_theoretical_solution_rin_rout(N,rin,rout)

f = 0.2;
fout = f;
w = 70;
p1 = rin/2/(1-f);
p2 = (1-f)/f*p1;
pout1 = rout/2/(1-fout);
pout2 = (1-fout)/fout*pout1;
%w = 70;
finh = 0.2;

A = (p1*(1-f) + (1-p2)*f)./ sqrt(f*(1-f)*(1-p1-p2).^2);
B = 1./sqrt(f*(1-f)*(1-p1-p2).^2);
C = ((1-f)*p1.*(1-p1)+ f*p2.*(1-p2))./(f*(1-f)*(1-p1-p2).^2);

Cout1 = sqrt(2)*erfinv(1-2*pout1);
Cout2 = sqrt(2)*erfinv(1-2*pout2);

Ninh = N*finh;
g=zeros(N,1);
g(1:Ninh) = -1;
g(Ninh+1:end) = 1;

E = @(x) (1+erf(x))/2;
F = @(x) exp(-x.^2)./pi^0.5+x.*(1+erf(x));
D = @(x) x.*F(x)+E(x);

options = optimset('Display','off','MaxIter',10^5,'MaxFunEvals',10000,'TolX',10^-8,'TolFun',10^-8);


S = @(x) [fout*F(x(2))-(1-fout)*F(x(1)); ...
    sum( (w*A.*g - B).*F(x(3).*(w*A.*g - B)/2)./(C + 4*(x(1)+x(2)).*( ...
    (fout*E(x(2))+(1-fout)*E(x(1)))./(fout * F(x(2))+(1-fout)*F(x(1))))./(Cout1 + Cout2)^2 ) );...
    sum( (2*(x(1) +x(2))^2 - (Cout1+Cout2)^2*C).*D(x(3)* (w*A.*g - B)/2)./(C + 4*(x(1)+x(2)).*( ...
    (fout*E(x(2))+(1-fout)*E(x(1)))./(fout * F(x(2))+(1-fout)*F(x(1))))./(Cout1 + Cout2)^2 ).^2 )];
j = 1;
while j < 5000
    j;
    x0 = [0 1 0]+0.1.*rand(1,3);
    [x,fval,exitflag] = fsolve(S, x0, options);
    if exitflag >0 && (x(1) + x(2))/(Cout1 + Cout2) >=0
        break;
    end
    j = j + 1;
end


if  exitflag >0
    fac = 8/(Cout1 + Cout2)^2*(fout*D(x(2))+(1-fout)*D(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)))^2;
    xx = 4*(x(1)+x(2))*((fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))))./(Cout1 + Cout2)^2;
    capacity = fac/N*sum( C.*D( x(3)* (w*A.*g - B)/2)./(C + xx ).^2 );
    Pcon = E((w*A.*g-B).* x(3)/2);
    Pconexc = Pcon(end);
    Pconinh = Pcon(1);
    kappa = 2*N*sum(C.*D((w*A.*g-B).* x(3)/2)./(C + xx ).^2)/(sum(A.*g.*F( x(3)* (w*A.*g - B)/2)./(C + xx)))^2;
    CV = sqrt(2*D((w*A.*g-B).* x(3)/2).*E((w*A.*g-B).* x(3)/2)./(F((w*A.*g-B).* x(3)/2)).^2 - 1);
    CVexc = CV(end);
    CVinh = CV(1);
    sumj= sqrt(2)*N*B./(xx + C)/sum(A.*g.*F(x(3).*(w.*A.*g - B)/2)./(xx+C));
    Jmean = g.*sumj.*F(x(3).*(w.*A.*g - B)/2)./E(x(3).*(w.*A.*g - B)/2)/sqrt(2);
    JmeanExc = Jmean(end);
    JmeanInh = Jmean(1);
    Wag = x(3).*(w.*A.*g - B)/2;
    JinhRange = (-1500 : 0.1 :0)/100;
    PropDensinh = w*exp(-(JinhRange*w/sqrt(2)/sumj + Wag(1)).^2)/sqrt(2*pi)/sumj/Pconinh;
    JexcRange = (0 : 0.1 :1500)/100;
    PropDensexc = w*exp(-(JexcRange*w/sqrt(2)/sumj - Wag(end)).^2) / sqrt(2*pi) /sumj / Pconexc;
%     figure
%     plot(JinhRange,PropDensinh,'color','b')
%     hold on
%     plot(JexcRange,PropDensexc,'color','r')
%     title('Weights distribution')
%     ylim([-0.005,0.4])
%     dim = [0.6 0.6 0.3 0.3];
%     %stt = {'plotting y=x^{var} ,',['with var =' num2str(var)]}
%     str = {['P_{con}^{inh}:',num2str(Pconinh)],['P_{con}^{exc}:',num2str(Pconexc)],['CV_{con}^{inh}:',num2str(CVinh)],['CV_{con}^{exc}:',num2str(CVexc)]};
%     annotation('textbox',dim,'String',str,'FitBoxToText','on');
%     axis square
else
    capacity = NaN;
    Pconinh = NaN;
    Pconexc = NaN;
end


