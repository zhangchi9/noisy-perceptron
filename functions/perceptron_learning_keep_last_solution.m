function [w,rout_perceptron] = perceptron_learning_keep_last_solution(X,Xp,f,rj,beta_post)

[N,m] = size(X);
Ninh = round(0.2*N);
beta = 0.1;
w1 = 70;
p1 = rj/2/(1-f);
p2 = rj/2/f;

X1 = (1 - p2*ones(1,m)).*X + (p1*ones(1,m)).*(1-X) ;
X2 = (p2 * ones(1,m)).*(1 - p2 * ones(1,m)).*X + (p1 * ones(1,m)).*(1 - p1 * ones(1,m)).*(1-X);

w= w1*[-ones(Ninh,1);ones(N-Ninh,1)];
g=[-ones(1,Ninh),ones(1,N-Ninh)]';

n_epochs = 50*10^5;
%rout_perceptron = nan(1,round(n_epochs/sample_interval));
%w_iteration = nan(N,round(n_epochs/sample_interval));
% Std_W = 1/N.*sqrt(N*int_noise^2 + sum(( W.^2 * ones(1,m)).*X2));
% Mean_W = Xp.*( W'*X1/N - 1);
% p1p2_W = 0.5*(1-erf(Mean_W./Std_W./sqrt(2)));
% rout_W = (1-f)*max(p1p2_W(Xp==-1)) + f*max(p1p2_W(Xp==1));
for epoch = 1: n_epochs
    for i = 1 : m
        x = X(:,i);
        errx = x;
        ind1 = find(x==0);
        ind1=ind1(randperm(length(ind1),round(length(ind1)*p1)));
        ind2 = find(x==1);
        ind2=ind2(randperm(length(ind2),round(length(ind2)*p2)));
        errx(ind1)=1;
        errx(ind2)=0;
        output = w'*errx/N - 1 + randn*beta_post/N^0.5 >= 0;
        if Xp(i)~=(2*output-1)
            w = w + beta.*Xp(i).*errx;
            w(w.*g<0)=0;
            delta_w = (sum(w.*g) - N*w1) /N;
            w = w - g.*delta_w;
            w(w.*g<0)=0;
        end
    end
    
    Std = 1/N.*sqrt(N*beta_post^2 + sum(( w.^2 * ones(1,m)).*X2));
    Mean = Xp.*( w'*X1/N - 1);
    p1p2 = 0.5*(1-erf(Mean./Std./sqrt(2)));
    rout_perceptron = (1-f)*max(p1p2(Xp==-1)) + f*max(p1p2(Xp==1));
end
end
