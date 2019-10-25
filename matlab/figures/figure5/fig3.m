clear
clc
close all
load dynamics_data_load_0.85_add_noise.mat

rin_range = rin_range_ind;

figure,imagesc(a_range,rin_range,CV_ISI)
hold on 
[C,h] = contour(a_range,rin_range,CV_ISI,[0.7 0.7], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
axis square
title('CV of ISI')
caxis([0.5,0.85])
colorbar 

figure,imagesc(a_range,rin_range,SPKS_COR)
hold on 
[C,h] = contour(a_range,rin_range,SPKS_COR,[0.15 0.15], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
axis square
title('SPKS COR')
colorbar 
caxis([0,0.2])

figure,imagesc(a_range,rin_range,COR_I)
hold on 
[C,h] = contour(a_range,rin_range,COR_I,[0.4 0.4], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
axis xy
xlabel('Postsynaptic noise strength, \beta_{post}')
ylabel('Spiking error probability, r ')
yticks(-12:2:-2)
yticklabels({'2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}'})
axis square
title('COR I')
colorbar 
caxis([0,1])
