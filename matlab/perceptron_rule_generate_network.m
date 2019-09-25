function perceptron_rule_generate_network(N,TrialNum)

a = 40;

rj = 2^-6;

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

rout = rj;
f = 0.2;
fout = f;

alpha = 0.1368;

p1 = rj/2/(1-f)*ones(N,1);
p2 = (1-f)*p1/f;
m = round(alpha*N);
X =rand(N,m+1)<f; %Input patterns
XX = X(:,1:end-1);

W = nan(N,N);
routmin = nan(1,N);

scratch_dir = ['/home/zhang.chi9/matlabtmp/', num2str(randi(10^9))];
mkdir(scratch_dir);
pc = parcluster('local');
pc.JobStorageLocation = scratch_dir;
parpool(pc,pc.NumWorkers);
parfor i=1:N
    i
    Xp = 2*X(i,2:end)-1;
    [W(:,i),routmin(i)] = perceptron_learning_with_noise(XX,Xp,f,rj,a,rout);
end
delete(gcp)
save(sprintf('/home/zhang.chi9/research/logscale/perceptron_capacity/N_%d_TrialNum_%d.mat',N,TrialNum))

end

function [w_return,routmin] = perceptron_learning_with_noise(X,Xp,f,rj,int_noise,rout)

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

n_epochs = 5*10^5;
routmin = 1000000;
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
        output = w'*errx/N - 1 + randn*int_noise/N^0.5 >= 0;
        if Xp(i)~=(2*output-1)
            w = w + beta.*Xp(i).*errx;
            w(w.*g<0)=0;
            delta_w = (sum(w.*g) - N*w1) /N;
            w = w - g.*delta_w;
            w(w.*g<0)=0;
        end
    end
    
    Std = 1/N.*sqrt(N*int_noise^2 + sum(( w.^2 * ones(1,m)).*X2));
    Mean = Xp.*( w'*X1/N - 1);
    p1p2 = 0.5*(1-erf(Mean./Std./sqrt(2)));
    rout_perceptron = (1-f)*max(p1p2(Xp==-1)) + f*max(p1p2(Xp==1));
    if rout_perceptron <= rout
        routmin = rout_perceptron;
        w_return = w;
        break
    elseif min(routmin,rout_perceptron ) < routmin
        routmin = rout_perceptron;
        w_return = w;
    end
    
end
% toc
%figure,hold on
%plot(sample_interval:sample_interval:n_epochs,rout_perceptron,'b')
%plot(sample_interval:sample_interval:n_epochs,rout_W*ones(1,length(rout_perceptron)),'r--')
%plot(sample_interval:sample_interval:n_epochs,rout*ones(1,length(rout_perceptron)))

% figure,plot(W,w,'o'),hold on
% plot(1.1*[min(W),max(W)],1.1*[min(W),max(W)])
% xlim(1.1*[min(W),max(W)])
% ylim(1.1*[min(W),max(W)])
% axis square
end
