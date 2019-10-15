
clear
clc
rin_ind_range =  -2:-0.5:-12;
a_range = [0.1,5:5:100];


for  i = 1: length(rin_ind_range)
    for j = 1:length(a_range)
        rin_ind = rin_ind_range(i);
        a = a_range(j);
    n_networks(i,j) = length(dir(['/scratch/zhang.chi9/perceptron/data/network_load_at_numerical_capacity/rin_',...
    num2str(rin_ind),'a_',num2str(a),'TrialNum_*.mat']));
    end
end

figure,imagesc(a_range,rin_ind_range,n_networks)
axis xy
title('new')