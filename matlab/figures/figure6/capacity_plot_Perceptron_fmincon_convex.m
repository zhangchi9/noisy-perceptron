% plot the capacity of perceptron, fmincon, convex
clear
clc
loadfac_range = 0.5:0.1:1.5;
perceptron_learning_rule = nan(1, length(loadfac_range));
fmincon_rule = nan(1, length(loadfac_range));
errbar_perceptron = nan(1, length(loadfac_range));
errbar_fmincon = nan(1, length(loadfac_range));
for i = 1 : length(loadfac_range)
    filename = ['/home/chi/Dropbox/Research/Perceptron_learning/data/perceptron_400/memory_load_',num2str(loadfac_range(i)),'.mat'];
    load(filename)
    perceptron_learning_rule(i) = sucprob_w;
    fmincon_rule(i) = sucprob_W;
    errbar_perceptron(i) = std(feasible_w(exitflag>0))/sqrt(length(feasible_w(exitflag>0)));
    errbar_fmincon(i) = std(feasible_W(exitflag>0))/sqrt(length(feasible_W(exitflag>0)));
end
load_fit = loadfac_range*0.1872;
figure,errorbar([0,loadfac_range*0.1872,0.4],[1,perceptron_learning_rule,0],[0,errbar_perceptron,0],'c'),hold on 
errorbar([0,loadfac_range*0.1872,0.4],[1,fmincon_rule,0],[0,errbar_fmincon,0],'m')
plot([0.1872,0.1872],[0,1],'k--')
legend('Perceptron','fmincon','replica')
axis square
    