function cover_theoretical()
N = 1000;
m_max=2.5*N;
C=zeros(1,m_max);
P=zeros(1,m_max);
if N<=50
    for m=1:m_max
        k_max=min([m-1,N-1]);
        C(m)=2.*sum(gamma(m-1+1)./gamma((0:k_max)+1)./gamma(m-1-(0:k_max)+1));
        P(m)=C(m)./2^m;
    end
else
    P=(1+erf((2*N-(1:m_max))./(2.*(1:m_max)).^0.5))./2;
end
plot((1:m_max)./N, P)
hold on 
plot((1:m_max)./N, P.^N)
ylim([0,1])
grid minor