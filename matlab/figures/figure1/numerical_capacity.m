clear
clc
close all
rin = -6;
a = 40;
loadfac = [1.5:-0.1:0.5];
prob_learning = nan(3,length(loadfac)+3,3);
N = [200, 400 800];
addpath('/home/chi/Dropbox/Research/Project_Perceptron/code/functions')
for i = 1 : length(N)
        for k = 1:length(loadfac)
            filename = ['/home/chi/Dropbox/Research/Project_Perceptron/data/figure1data/N_',num2str(N(i)),'_rj_',num2str(rin),'_a_',num2str(a),'_loadfac_',num2str(loadfac(k)),'.mat'];
            if isfile(filename)
                N(i)
                file_variable = load(filename);
                epsilonsum = file_variable.epsilonsum;
                exitflag = file_variable.exitflag;
                prob_learning(i,k+1,1) = file_variable.m/N(i);
                prob_learning(i,k+1,2) = mean(epsilonsum(exitflag>0) < 10^-4);
                prob_learning(i,k+1,3) = std(epsilonsum(exitflag>0) < 10^-4)/sqrt(numel(epsilonsum(exitflag>0)));
            end
        end
end
c = theoretical_solution(a,0,2.^rin,0.2);
figure,
prob_learning(:,length(loadfac)+2,1) = 0.07;
prob_learning(:,length(loadfac)+2,2) = 1;
prob_learning(:,length(loadfac)+2,3) = 0;
prob_learning(:,length(loadfac)+3,[1,3]) = 0;
prob_learning(:,length(loadfac)+3,2) = 1;
prob_learning(:,1,1) = 0.4;
prob_learning(:,1,[2,3]) = 0;

errorbar(prob_learning(1,:,1),prob_learning(1,:,2),prob_learning(1,:,3))
hold on 
errorbar(prob_learning(2,:,1),prob_learning(2,:,2),prob_learning(2,:,3))
errorbar(prob_learning(3,:,1),prob_learning(3,:,2),prob_learning(3,:,3))
plot([c,c],[0,1],'--','color','k')
xlim([0,0.4])
axis square
legend('N=200','N=400','N=800','N->\infty')

get_numerical_capacity([prob_learning(1,:,1);prob_learning(1,:,2)]')
get_numerical_capacity([prob_learning(2,:,1);prob_learning(2,:,2)]')
get_numerical_capacity([prob_learning(3,:,1);prob_learning(3,:,2)]')
