clear
clc
close all
load retrieval_prob_length_map_load85_same_noise.mat
load 'D:\Dropbox\research\Project_Perceptron\codes\figures\figure2\fig2data.mat'
rin_range_ind =  -2:-0.5:-12;
a_range = [0.1,5:5:100];
noise_added = 0.1;%[0.1,2,4,5:5:150];
N = 400;
retrieval_length_map = nan(length(rin_range_ind),length(a_range));
retrievl_prob_map = nan(length(rin_range_ind),length(a_range));
for i = 1:length(rin_range_ind)
    for j = 1:length(a_range)
        index = 1;%find(noise_added == a_range(j));
        retrieval_length_map(i,j) = retrieval_length(i,j,index);
        retrievl_prob_map(i,j) = retrieval_prob(i,j,index);
    end
end
capacity = capacity([end:-1:1],:);
I = I([end:-1:1],:);

figure,imagesc(a_range,rin_range_ind,retrieval_length_map)
axis xy
colorbar
title('retrieval length')
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
axis square

figure,imagesc(a_range,rin_range_ind,retrievl_prob_map)
axis xy
colorbar
title('retrieval probability')
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
axis square

% figure,imagesc(a_range,rin_range_ind,retrieval_length_map/N./capacity/0.85)
% axis xy
% colorbar
% title('retrieval fraction')
% xlabel('Postsynaptic noise strength, \beta_{post}')
% ylabel('Spiking error probability, r ')
% yticks(-12:2:-2)
% yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
% axis square
% 
% figure,imagesc(a_range,rin_range_ind,retrieval_length_map.*I./capacity)
% axis xy
% colorbar
% title('retrievl Information')
% xlabel('Postsynaptic noise strength, \beta_{post}')
% ylabel('Spiking error probability, r ')
% yticks(-12:2:-2)
% yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
% caxis([0,25])
% axis square