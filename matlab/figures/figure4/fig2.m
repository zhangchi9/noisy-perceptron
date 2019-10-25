
clear
clc
close all
beta = 0;
a_range = [0.1,5:5:100];
%rin_range = [0.01:0.01:0.2];
rin_range = 2.^([-2:-0.5:-12]);
f = 0.2;
model = 'homo';

for i = 1 : length(rin_range)
    for j = 1 : length(a_range)
        [capacity(i,j),exitflag,Pconinh(i,j),Pconexc(i,j),CVinh(i,j),CVexc(i,j)] = theoretical_solution(a_range(j),beta,rin_range(i),f);
    end
end

fout = 0.2;

Is = fout*(rin_range/2/fout.*log(rin_range/2/fout) + (1-rin_range/2/fout).*log(1-rin_range/2/fout) -log(fout)) + ...
    (1-fout)*(rin_range/2/(1-fout).*log(rin_range/2/(1-fout)) + (1-rin_range/2/(1-fout)).*log(1-rin_range/2/(1-fout)) -log((1-fout)));
I = capacity.* (Is'*ones(1,length(a_range)));

rin_ind_range = (log2(rin_range));

figure,imagesc(a_range,rin_ind_range,I)
hold on 
[C,h] = contour(a_range,rin_ind_range,I,[0.05 0.1 0.2], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('Information')
axis square
caxis([0.1 0.2])


load('../figure4/retrieval_prob_length_map_load85_no_noise.mat')
figure,imagesc(a_range,rin_ind_range,retrieval_prob.*I)
hold on 
[C,h] = contour(a_range,rin_ind_range,retrieval_prob.*I,[0.05 0.1], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('Information')
axis square
caxis([0 0.15])

figure,imagesc(a_range,rin_ind_range,capacity)
hold on 
[C,h] = contour(a_range,rin_ind_range,capacity,[0.3 0.3], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
[C,h] = contour(a_range,rin_ind_range,capacity,[0.2 0.2], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('Capacity')
axis square
caxis([0.1 0.4])

figure,imagesc(a_range,rin_ind_range,Pconinh)
hold on 
[C,h] = contour(a_range,rin_ind_range,Pconinh,[0.25,0.56], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('Pconinh')
axis square
caxis([0 0.7])

figure,imagesc(a_range,rin_ind_range,Pconexc)
hold on 
[C,h] = contour(a_range,rin_ind_range,Pconexc,[0.1,0.19], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('Pconexc')
axis square
caxis([0 0.3])

figure,imagesc(a_range,rin_ind_range,CVinh)
hold on 
[C,h] = contour(a_range,rin_ind_range,CVinh,[0.78,0.96], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('CVinh')
axis square
caxis([0.7 0.9])

figure,imagesc(a_range,rin_ind_range,CVexc)
hold on 
[C,h] = contour(a_range,rin_ind_range,CVexc,[0.85,1.1], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
title('CVexc')
axis square
caxis([0.7 0.9])

