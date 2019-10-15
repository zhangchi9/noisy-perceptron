close all
clear
load 'single_retriveal_length_prob,mat.mat'
figure
xlabel('noise added in retrieval/noise added in training')
yyaxis left
errorbar((5:5:100)/40,prob,prob_std)
ylabel('retrieval prob')
ylim([0,1])
yyaxis right
errorbar((5:5:100)/40,retriveal_length/max(retriveal_length),retriveal_length_std/max(retriveal_length))
ylabel('normalized retrieval length')
axis square
ylim([0,1])
