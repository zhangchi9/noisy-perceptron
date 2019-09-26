% This function will calculate the robust model theoretical_solution,
% f: presynaptic neuron firing rate, p presynaptic neuron error rate
% fout: postsynaptic neuron firing rate, pout presynaptic neuron error rate
% w: L1 norm, finh: inhibitory neuron percentage
% All the probability are independent of j

function [capacity,exitflag,Pconinh,Pconexc,CVinh,CVexc] = theoretical_solution_all_cases(error_in,error_out,model)

N = 1000;
f = 0.2;
fout = f;
w = 70;
finh = 0.2;
Ninh = N*finh;

pout1 = error_out/2/(1-fout);
pout2 = (1-fout)/fout*pout1;
Cout1 = sqrt(2)*erfinv(1-2*pout1);
Cout2 = sqrt(2)*erfinv(1-2*pout2);

E = @(x) (1+erf(x))/2;
F = @(x) exp(-x.^2)./pi^0.5+x.*(1+erf(x));
D = @(x) x.*F(x)+E(x);

g=zeros(N,1);
g(1:Ninh) = -1;
g(Ninh+1:end) = 1;

switch model
    case 1
        options = optimset('Display','off','MaxIter',10^5,'MaxFunEvals',10000,'TolX',10^-8,'TolFun',10^-8);
        
        S = @(x) [fout*F(x(2))-(1-fout)*F(x(1)); ...
            finh*(1+w*f)*F((2*x(4)*f-x(3))/2/sqrt(f*(1-f)))+ (1-finh)*(1-w*f)*F((-2*x(4)*f-x(3))/2/sqrt(f*(1-f)));...
            -finh*f/sqrt(f*(1-f))*F((2*x(4)*f-x(3))/2/sqrt(f*(1-f)))+ (1-finh)*f/sqrt(f*(1-f))*F((-2*x(4)*f-x(3))/2/sqrt(f*(1-f))) ...
            - 8*(x(1)+x(2))/(Cout1+Cout2)^2/error_in^2 * (0.5*w*x(3)+x(4))* (fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)));...
            finh*D((2*x(4)*f-x(3))/2/sqrt(f*(1-f)))+(1-finh)*D((-2*x(4)*f-x(3))/2/sqrt(f*(1-f)))-16/(Cout1+Cout2)^2/error_in^2 * (0.5*w*x(3)+x(4))^2*...
            ((fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))))^2];
        
        j = 1;
        while j < 5000
            j
            x0 = [0 1 0 0]+0.1.*rand(1,4);
            [x,~,exitflag] = fsolve(S, x0, options);
            if exitflag >0 && (x(1) + x(2))* (0.5*w*x(3)+x(4)) > 0 &&(x(1) + x(2))/(Cout1 + Cout2) >=0
                break;
            end
            j = j + 1;
        end
        if  exitflag >0
            fac = 16/(Cout1 + Cout2)^2/error_in^2*(fout*D(x(2))+(1-fout)*D(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)))^2;
            capacity = fac* (0.5*w*x(3)+x(4))^2;
            Pconexc = E(-(x(3)+2*x(4)*f)/2/sqrt(f*(1-f)));
            Pconinh = E(-(x(3)-2*x(4)*f)/2/sqrt(f*(1-f)));
            CV = g.*sqrt(2*D(-(x(3)+2*x(4)*f*g)/2/sqrt(f*(1-f))).*E(-(x(3)+2*x(4)*f*g)/2/sqrt(f*(1-f)))...
                ./(F(-(x(3)+2*x(4)*f*g)/2/sqrt(f*(1-f)))).^2 - 1);
            CVinh = -CV(1);
            CVexc = CV(end);
        end
        
    case 2
        options = optimset('Display','off','MaxIter',10^5,'MaxFunEvals',10000,'TolX',10^-8,'TolFun',10^-8);
        
        S = @(x) [fout*F(x(2))-(1-fout)*F(x(1)); ...
            finh*f/sqrt(f*(1-f))*F(-(x(3) - x(4)*f*(1-1/f/w))/sqrt(f*(1-f)))...
            - 8*(x(1)+x(2))/(Cout1+Cout2)^2/error_in* (0.5*w*x(3)+x(4))*(1-1/f/w)*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)));...
            (1-finh)*f/sqrt(f*(1-f))*F(-(x(3) + x(4)*f*(1+1/f/w))/sqrt(f*(1-f)))...
            - 8*(x(1)+x(2))/(Cout1+Cout2)^2/error_in* (0.5*w*x(3)+x(4))*(1+1/f/w)*(fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)));
            finh*D(-(x(3) - x(4)*f*(1-1/f/w))/sqrt(f*(1-f))) + (1-finh)*D(-(x(3) + x(4)*f*(1+1/f/w))/sqrt(f*(1-f)))...
            - 64/(Cout1+Cout2)^2/w/f/error_in* (0.5*w*x(3)+x(4))^2*((fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))))^2];
        
        j = 1;
        while j < 5000
            j
            x0 = [0 1 0 0]+0.1.*rand(1,4);
            [x,~,exitflag] = fsolve(S, x0, options);
            if exitflag >0 && (x(1) + x(2))* (0.5*w*x(3)+x(4)) > 0 &&(x(1) + x(2))/(Cout1 + Cout2) >=0
                break;
            end
            j = j + 1;
        end
        
        if  exitflag >0
            fac = 64/(Cout1+Cout2)^2/w/f/error_in*(fout*D(x(2))+(1-fout)*D(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)))^2;
            capacity = fac* (0.5*w*x(3)+x(4))^2;
            Pconexc = E(-(x(3)+x(4)/w+x(4)*f)/sqrt(f*(1-f)));
            Pconinh = E(-(x(3)+x(4)/w-x(4)*f)/sqrt(f*(1-f)));
            CV = g.*sqrt(2*D(-(x(3)+x(4)/w+x(4)*f*g)/sqrt(f*(1-f))).*E(-(x(3)+x(4)/w+x(4)*f*g)/sqrt(f*(1-f)))...
                ./(F(-(x(3)+x(4)/w+x(4)*f*g)/sqrt(f*(1-f)))).^2 - 1);
            CVinh = -CV(1);
            CVexc = CV(end);
        end
        
    case 3
        
        p1 = error_in/2/(1-f);
        p2 = (1-f)/f*p1;
        
        Aj = p1*(1-f) + (1-p2)*f;
        Cj = (1-f)*p1.*(1-p1)+ f*p2.*(1-p2);
        Dj = f*(1-f)*(1-p1-p2).^2;
       
        options = optimset('Display','off','MaxIter',10^5,'MaxFunEvals',10000,'TolX',10^-8,'TolFun',10^-8);
        
        S = @(x) [fout*F(x(2))-(1-fout)*F(x(1)); ...
            sum( (w*Aj.*g - 1)*sqrt(Dj).*F(x(3).*(w*Aj.*g - 1)/2/sqrt(Dj))./(Cj + Dj*4*(x(1)+x(2)).*( ...
            (fout*E(x(2))+(1-fout)*E(x(1)))./(fout * F(x(2))+(1-fout)*F(x(1))))./(Cout1 + Cout2)^2 ) );...
            sum( (Dj*(x(1) +x(2))^2/(Cout1+Cout2)^2 - Cj/2).*D(x(3).*(w*Aj.*g - 1)/2/sqrt(Dj))*Dj./(Cj + Dj*4*(x(1)+x(2)).*( ...
            (fout*E(x(2))+(1-fout)*E(x(1)))./(fout * F(x(2))+(1-fout)*F(x(1))))./(Cout1 + Cout2)^2 ).^2 )];
        
        j = 1;
        while j < 500
            j
            x0 = [0 1 0]+0.1.*rand(1,3);
            [x,~,exitflag] = fsolve(S, x0, options);
            tmp = sum( Aj.*g * sqrt(Dj).*F(x(3).*(w*Aj.*g - 1)/2/sqrt(Dj))./(Cj + Dj*4*(x(1)+x(2)).*( ...
                (fout*E(x(2))+(1-fout)*E(x(1)))./(fout * F(x(2))+(1-fout)*F(x(1))))./(Cout1 + Cout2)^2 ) );
            if exitflag >0 && (x(1) + x(2))*tmp > 0 &&(x(1) + x(2))/(Cout1 + Cout2) >=0
                break;
            end
            j = j + 1;
        end
        
        if  exitflag >0
            fac = 8/(Cout1 + Cout2)^2*(fout*D(x(2))+(1-fout)*D(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1)))^2;
            xx = 4*(x(1)+x(2))*((fout*E(x(2))+(1-fout)*E(x(1)))/(fout * F(x(2))+(1-fout)*F(x(1))))./(Cout1 + Cout2)^2;
            capacity = fac/N*sum( Cj*Dj.*D( x(3)* (w*Aj.*g - 1)/2/sqrt(Dj))./(Cj + Dj*xx ).^2 );
            Pcon = E((w*Aj.*g-1).* x(3)/2/sqrt(Dj));
            Pconinh = Pcon(1);
            Pconexc = Pcon(end);
            
            CV = g.*sqrt(2*D((w*Aj.*g-1).* x(3)/2/sqrt(Dj)).*E((w*Aj.*g-1).* x(3)/2/sqrt(Dj))...
                ./(F((w*Aj.*g-1).* x(3)/2/sqrt(Dj))).^2 - 1);
            CVinh = -CV(1);
            CVexc = CV(end);
            
            %kappa = 2*N*sum(Cj.*D((w*Aj.*g-B).* x(3)/2)./(Cj + xx ).^2)/(sum(Aj.*g.*F( x(3)* (w*Aj.*g - B)/2)./(Cj + xx)))^2;
            %CV = sqrt(2*D((w*Aj.*g-B).* x(3)/2).*E((w*Aj.*g-B).* x(3)/2)./(F((w*Aj.*g-B).* x(3)/2)).^2 - 1);
            %sumj= sqrt(2)*N*B./(xx + Cj)/sum(Aj.*g.*F(x(3).*(w.*Aj.*g - B)/2)./(xx+Cj));
            %Jmean = g.*sumj.*F(x(3).*(w.*Aj.*g - B)/2)./E(x(3).*(w.*Aj.*g - B)/2)/sqrt(2);
            
            %Wag = g.*x(3).*(w.*Aj.*g - B)/2;
            % a is Jrange,b is sumj, c is pcon, d is Wag
            %PropDensinh = @(a,b,c,d) w*exp(-(a*w/sqrt(2)/b - d).^2)/sqrt(2*pi)/b/c;
            %     figure,
            %     subplot(4,2,1),plot(rin(1:Ninh),Pcon(1:Ninh))
            %     [~,pos]=min(abs(Pcon(1:Ninh) - mean(Pcon(1:Ninh))));
            %     hold on, plot(rin(pos),mean(Pcon(1:Ninh)),'o')
            %
            %     subplot(4,2,2),plot(rin((Ninh+1):N),Pcon((Ninh+1):end))
            %     [~,pos]= min(abs(Pcon((Ninh+1):N) - mean(Pcon((Ninh+1):N))));
            %     hold on, plot(rin(pos+Ninh),mean(Pcon((Ninh+1):N)),'o')
            %
            %     subplot(4,2,3),plot(rin(1:Ninh),CV(1:Ninh))
            %     [~,pos]=min(abs(CV(1:Ninh) - mean(CV(1:Ninh))));
            %     hold on, plot(rin(pos),mean(CV(1:Ninh)),'o')
            %
            %     subplot(4,2,4),plot(rin((Ninh+1):N),CV((Ninh+1):end))
            %     [~,pos]=min(abs(CV((Ninh+1):N) - mean(CV((Ninh+1):N))));
            %     hold on, plot(rin(pos+Ninh),mean(CV((Ninh+1):N)),'o')
            %
            %     subplot(4,2,5),plot(rin(1:Ninh),-Jmean(1:Ninh))
            %     subplot(4,2,6),plot(rin((Ninh+1):N),Jmean((Ninh+1):end))
            %     subplot(4,2,7),plot(rin(1:Ninh),-Jmean(1:Ninh).*Pcon(1:Ninh))
            %     subplot(4,2,8),plot(rin((Ninh+1):N),Jmean((Ninh+1):end).*Pcon((Ninh+1):end))
            %
        else
            capacity = NaN;
            Pconinh = NaN;
            Pconexc = NaN;
        end
    otherwise
        disp('other value')
end





