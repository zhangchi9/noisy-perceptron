clear
clc
close all
rin = -6;
a = 40;
loadfac = [1.5:-0.1:0.5];
prob_learning = nan(3,length(loadfac)+3,3);
N = 800;

for k = 1:length(loadfac)
    filename = ['N_',num2str(N),'_rj_',num2str(rin),'_a_',num2str(a),'_loadfac_',num2str(loadfac(k)),'.mat'];
    if isfile(filename)
        file_variable = load(filename);
        
        prob_learning(1,k+1,1) = file_variable.m/N;
        prob_learning(1,k+1,2) = file_variable.prob_learning;
        tmp = file_variable.epsilonsum;
        prob_learning(1,k+1,3) = std(tmp(:)<0.001)/sqrt(numel(tmp(:)));
    end
end

for k = 1:length(loadfac)
    filename = ['N_',num2str(N),'_rj_',num2str(rin),'_beta_',num2str(a/14),'_loadfac_',num2str(loadfac(k)),'.mat'];
    if isfile(filename)
        file_variable = load(filename);
        
        prob_learning(2,k+1,1) = file_variable.m/N;
        prob_learning(2,k+1,2) = file_variable.prob_learning;
        tmp = file_variable.epsilonsum;
        prob_learning(2,k+1,3) = std(tmp(:)<0.001)/sqrt(numel(tmp(:)));
    end
end

c = theoretical_solution(a,0,2.^rin,0.2);
figure,
errorbar(prob_learning(1,:,1),prob_learning(1,:,2),prob_learning(1,:,3))
hold on
errorbar(prob_learning(2,:,1),prob_learning(2,:,2),prob_learning(2,:,3))
plot([c,c],[0,1],'--','color','k')
xlim([0,0.4])
axis square
legend('\beta_{int} = 40','\beta_{syn} = \beta_{int}^2/wf','N->\infty')
