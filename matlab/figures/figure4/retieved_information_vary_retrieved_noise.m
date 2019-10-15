clear
clc

load('/home/chi/Dropbox/Research/Project_Perceptron/results/retrieval_prob_length_map_load_at_numerical_capacity_varied_retried_noise_fixed_rj.mat')
load('/home/chi/Dropbox/Research/Project_Perceptron/results/m_capacity.mat')

capacity = m_capacity(9,:)/400;
fout = 0.2;

rin_range = 2.^rin_range_ind;

Is = fout*(rin_range/2/fout.*log(rin_range/2/fout) + (1-rin_range/2/fout).*log(1-rin_range/2/fout) -log(fout)) + ...
    (1-fout)*(rin_range/2/(1-fout).*log(rin_range/2/(1-fout)) + (1-rin_range/2/(1-fout)).*log(1-rin_range/2/(1-fout)) -log((1-fout)));
I = ones(length(a_range),1)*(capacity.* (Is'*ones(1,length(a_range))));

rin_ind_range = (log2(rin_range));

retrieval_prob = squeeze(retrieval_prob)';
retrieval_length = squeeze(retrieval_length)'/400;

figure,imagesc(a_range,a_range,retrieval_prob)
hold on 
[C,h] = contour(a_range,a_range,retrieval_prob,[0.5 0.5], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post} learning')
ylabel('Postsynaptic noise strength, \beta_{post} retrieval')
title('retrieval probability')
axis square
caxis([0 1])

figure,imagesc(a_range,a_range,retrieval_length./capacity)
hold on 
[C,h] = contour(a_range,a_range,retrieval_length./capacity,[0.5 0.5], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post} learning')
ylabel('Postsynaptic noise strength, \beta_{post} retrieval')
title('retrieval length fraction')
axis square
caxis([0 1])

figure,imagesc(a_range,a_range,retrieval_prob.*I)
hold on 
[C,h] = contour(a_range,a_range,retrieval_prob.*I,[0.04 0.04], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post} learning')
ylabel('Postsynaptic noise strength, \beta_{post} retrieval')
title('retrieval probability * information')
axis square
caxis([0 0.05])

figure,imagesc(a_range,a_range,retrieval_length./capacity.*I)
hold on 
[C,h] = contour(a_range,a_range,retrieval_length./capacity.*I,[0.04 0.04], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post} learning')
ylabel('Postsynaptic noise strength, \beta_{post} retrieval')
title('retrieval information')
axis square
caxis([0 0.05])

