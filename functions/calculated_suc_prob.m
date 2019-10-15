function calculated_suc_prob(rj,beta_post,loadfac)

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

N=100;
rout = rj; 
f = 0.2;
fout = f;
totol_training = 200;

alpha = theoretical_solution(beta_post,0,rj,f);
m = round(loadfac* alpha * N);
p1 = rj/2/(1-f)*ones(N,1);
p2 = (1-f)*p1/f;

W = nan(N,totol_training);
epsilonsum = nan(1,totol_training);
exitflag = nan(1,totol_training);

filename = (['N_',num2str(N),'_rj_',num2str(rj),'_a_',num2str(beta_post),'_loadfac_',num2str(loadfac),'.mat']);

% scratch_dir = ['/home/zhang.chi9/matlabtmp/', num2str(randi(10^9))];
% mkdir(scratch_dir);
% pc = parcluster('local');
% pc.JobStorageLocation = scratch_dir;
% parpool(pc,pc.NumWorkers);
% par
for i = 1 : totol_training
    i
    X = rand(N,m)<f;
    Xp = 2*(rand(1,m)<fout)-1;
    
    [W(:,i),epsilonsum(i),exitflag(i)] = numerical_solution(X,Xp,f,rj,beta_post,rout,'fmincon');%find_W_single_trial(X,X1,X2,Xp,C,bate_int);
end
%delete(gcp)
prob_learning = mean(epsilonsum<10^-6);
save(filename)
end


