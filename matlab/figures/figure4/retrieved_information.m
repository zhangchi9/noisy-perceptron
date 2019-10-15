clear
clc

load('/home/chi/Dropbox/Research/Project_Perceptron/results/retrieval_prob_length_map_load_at_numerical_capacity_noise_all.mat')
%load('retrieval_prob_length_map_load85_no_noise.mat')
load('/home/chi/Dropbox/Research/Project_Perceptron/results/m_capacity.mat')

capacity = m_capacity/400;
fout = 0.2;

rin_range = 2.^rin_range_ind;

retrieval_prob = retrieval_prob(:,:,8);
retrieval_length = retrieval_length(:,:,8);

Is = fout*(rin_range/2/fout.*log(rin_range/2/fout) + (1-rin_range/2/fout).*log(1-rin_range/2/fout) -log(fout)) + ...
    (1-fout)*(rin_range/2/(1-fout).*log(rin_range/2/(1-fout)) + (1-rin_range/2/(1-fout)).*log(1-rin_range/2/(1-fout)) -log((1-fout)));
I = capacity.* (Is'*ones(1,length(a_range)));

rin_ind_range = (log2(rin_range));

figure,imagesc(a_range,rin_ind_range,retrieval_prob)
hold on 
[C,h] = contour(a_range,rin_ind_range,retrieval_prob,[0.5 0.5], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('retrieval probability')
axis square
caxis([0 1])

figure,imagesc(a_range,rin_ind_range,retrieval_length./m_capacity)
hold on 
[C,h] = contour(a_range,rin_ind_range,retrieval_length./m_capacity,[0.5 0.5], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('retrieval length fraction')
axis square
caxis([0 1])

figure,imagesc(a_range,rin_ind_range,retrieval_prob.*I)
hold on 
[C,h] = contour(a_range,rin_ind_range,retrieval_prob.*I,[0.05 0.05], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('Information in completely retrieved squences per N^2')
axis square
caxis([0 0.07])


V = retrieval_length./m_capacity.*I;
% K = ones(3)/9; 
% V = conv2(V,K,'same');

figure,imagesc(a_range,rin_ind_range,V)
hold on 
[C,h] = contour(a_range,rin_ind_range,V,[0.05 0.05], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('Retrieved sequence information per N^2')
axis square
caxis([0 0.07])


