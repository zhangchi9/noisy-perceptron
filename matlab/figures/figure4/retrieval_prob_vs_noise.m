clear
clc
load('/home/chi/Dropbox/Research/Project_Perceptron/results/retrieval_prob_length_map_load_at_numerical_capacity_noise_all.mat')
%load('retrieval_prob_length_map_load85_no_noise.mat')
load('/home/chi/Dropbox/Research/Project_Perceptron/results/m_capacity.mat')

capacity = m_capacity/400;
fout = 0.2;

rin_range = 2.^rin_range_ind;

retrieval_prob = retrieval_prob(9,:,8);
retrieval_length = retrieval_length(9,:,8); 

figure,yyaxis left, plot(a_range,retrieval_prob)
title('retrieval prob vs. noise')
ylabel('retrieval prob')
axis square

hold on,yyaxis right,plot(a_range,retrieval_length./m_capacity(9,:))
title('retrieval length vs. noise')
ylabel('retrieval length fraction')
xlabel('noise added in learning')
axis square