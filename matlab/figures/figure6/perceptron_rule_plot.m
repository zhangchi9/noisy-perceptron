% this script will plot the perceptron learning rule process

clear
clc
load('/home/chi/Dropbox/Perceptron_learning/data/TrialNum_106.mat')
feasible = epsilonsum<10^-4;
learned = routmin < rj;
ind = find(feasible == 1 & learned == 1);
W = W(:,ind(5));
w = w(:,ind(5));

Xp = 2*X(ind(5),2:end)-1;
X = XX;

Ninh = round(0.2*N);
w1 = 70;
f = 0.2;
fout = 0.2;
beta = 0.1;

noise_int = 40;
rj = 2^-6;
rout = 2^-6;
p1 = rj/2/(1-f);
p2 = rj/2/f;
pout1 = rout/2/(1-f);
pout2 = rout/2/f;


% feasible = 1;
% while feasible == 1
%     [W,epsilonsum,exitflag] = fmincon_with_noise(X,Xp,noise_int,rj,rout);
%     feasible = epsilonsum < 10^-4;
% end

%load test

X1 = (1 - p2*ones(1,m)).*X + (p1*ones(1,m)).*(1-X) ;
X2 = (p2 * ones(1,m)).*(1 - p2 * ones(1,m)).*X + (p1 * ones(1,m)).*(1 - p1 * ones(1,m)).*(1-X);

J= w1*[-ones(Ninh,1);ones(N-Ninh,1)];
g=[-ones(1,Ninh),ones(1,N-Ninh)]';

n_epochs = 100*10^4;
sample_interval = 1000;

rout_perceptron = nan(1,round(n_epochs/sample_interval));

% Std_W = 1/N.*sqrt(N*noise_int^2 + sum(( W.^2 * ones(1,m)).*X2));
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
        intrinsic_niose = randn*noise_int/N^0.5;
        output = J'*errx/N - 1 + intrinsic_niose >= 0;
        %output = Xp(i).*(w'*x/N-1)<= kappa/N^0.5;
        % flip a coin to decide if update the weiths
        if Xp(i)==-1 && output == 1 %&& rand < 1 - pout2
            J = J + beta.*Xp(i).*errx;
            J(J.*g<0)=0;
            delta_w = (sum(J.*g) - N*w1) /N;
            J = J - g.*delta_w;
            J(J.*g<0)=0;
            %         elseif Xp(i)== 1 && output == 1 && rand < pout1
            %             w = w + beta.*Xp(i).*x;
            %             w(w.*g<0)=0;
            %             delta_w = (sum(w.*g) - N*w1) /N;
            %             w = w - g.*delta_w;
            %             w(w.*g<0)=0;
            %         elseif Xp(i)== -1 && output == 0 && rand < pout2
            %             w = w + beta.*Xp(i).*x;
            %             w(w.*g<0)=0;
            %             delta_w = (sum(w.*g) - N*w1) /N;
            %             w = w - g.*delta_w;
            %             w(w.*g<0)=0;
        elseif Xp(i)==1 && output == 0 %&& rand < 1 - pout1
            J = J + beta.*Xp(i).*errx;
            J(J.*g<0)=0;
            delta_w = (sum(J.*g) - N*w1) /N;
            J = J - g.*delta_w;
            J(J.*g<0)=0;
        end
    end
    
    if mod(epoch,sample_interval) == 0
        epoch
        
        Std = 1/N.*sqrt(N*noise_int^2 + sum(( J.^2 * ones(1,m)).*X2));
        Mean = Xp.*( J'*X1/N - 1);
        p1p2 = 0.5*(1-erf(Mean./Std./sqrt(2)));
        rout_perceptron(round(epoch/sample_interval)) = (1-f)*max(p1p2(Xp==-1)) + f*max(p1p2(Xp==1));
        %         plot(W,w,'o'),hold on
        %         plot(1.1*[min(W),max(W)],1.1*[min(W),max(W)])
        %         xlim(1.1*[min(W),max(W)])
        %         text(-400,100,['iteration: ',num2str(epoch)])
        %         pause(0.1)
        %         clf
    end
end
%toc
figure,hold on
plot(sample_interval:sample_interval:n_epochs,rout_perceptron,'b')
%plot(sample_interval:sample_interval:n_epochs,rout_W*ones(1,length(rout_perceptron)),'r--')
plot(sample_interval:sample_interval:n_epochs,rout*ones(1,length(rout_perceptron)))

save('feasible_solved_by_perceptron.mat')




